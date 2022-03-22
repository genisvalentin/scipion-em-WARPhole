# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Genis Valentin Gese (genis.valentin.gese@ki.se)
# *
# * Karolinska Institutet
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'genis.valentin.gese@ki.se'
# *
# **************************************************************************
'''
Something to fix!
FATAL ERROR: Object 5750.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 5750.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 5750.outputMovies1 has no sampling rate!!!

Something to fix!
Check for new files after the correct fileTimeout
'''

"""
This protocol imports data pre-processed in WARP.
It will produce a set or particles, coordinates, micrographs, ctfmodels and movies.
This protocol can work in streaming mode.
"""

from pyworkflow.protocol import Protocol, params, Integer
from pwem.protocols.protocol_import.images import ProtImportImages
from pwem.objects import Micrograph, MovieAlignment, Movie, Particle
from pwem.objects.data import SetOfMicrographs, SetOfParticles, SetOfMovies, SetOfCoordinates, SetOfCTF
from pyworkflow.utils import Message
import pyworkflow.utils as pwutils
from .WARPimporter import WARPimporter
import time
import os
from pwem.protocols import EMProtocol
import pwem.emlib.metadata as md
import pyworkflow.protocol.constants as cons
import pickle
import itertools

class WARPholeImportParticles(EMProtocol):
    """
    This protocol imports cryo-EM data pre-processed by WARP from the goodparticles star file.
    """
    _label = 'Import WARP particles'
    _outputClassName = 'SetOfParticles'
    _importedMicrographsSet = list()
    _closed = False
    _readyToClose = False
    _stepsCheckSecsCounter = 4
    _lastImportedId = 0

    # -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('starFile', params.FileParam,
                      label='Star file',
                      important=True,
                      help="Select a *_data.star file from a\n"
                           "previous Relion execution."
                           "To detect if the input particles contains alignment "
                           "information, it is required to have the "
                           "optimiser.star file corresponding to the data.star")

        form.addParam('fileTimeout', params.FloatParam,
                      label='Load new particles after (sec): ',
                      default=0,
                      help="The star file is read every these many seconds, and new particles are imported. If no new particles are found, the import is considered to be finished. Set to zero for no streaming.")

        form.addParam('copyBinaries', params.BooleanParam,
                      default=False,
                      label='Copy binary files?',
                      important=False,
                      help="If no, the plugin will create symlinks to the imported binary files. If yes, the binary files will be copied into the scipion directory.")

        form.addParam('dosePerFrame', params.FloatParam,
                      label='Dose per frame',
                      default=0,
                      help="Set to zero if unknown")

        form.addParam('magnification', params.FloatParam,
                      label='Magnification',
                      default=0,
                      help="Set to zero if unknown")

        form.addParam('moviePixelSize', params.FloatParam,
                      label='Movie pixel size (unbinned)',
                      help="Pixel size of the unbinned movies. This is needed in case some binning is done in WARP.",
                      default=1)

        form.addParam('doImportAlignedMovies',
                      params.BooleanParam,
                      label='Wait for movies alignment star files?',
                      default=False,
                      help="If yes, the plugin will wait for movie alignments to be exported by WARP, and output a set of aligned micrographs. These are produced by seleting generate files for relion polishing in the export micrograph list menu of WARP. Useful to run relion bayesian polishing.")

        form.addParam('movieTimeout', params.FloatParam,
                      label='Wait for movie alignment files after import? (sec)',
                      condition='(doImportAlignedMovies == True)',
                      default=72000,
                      help="After the particle import is finished, what these many seconds for the movie motion correction star files to be available. Set to zero if the files are already available.")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        #Create empty sets of particles, coordinates, micrographs, ctfs and movies.
        self._insertFunctionStep('prepareImporterStep')
        #Prepare a step for closing the particle sets
        self.closeSetsId = self._insertFunctionStep('closeSetsStep', wait=True)
        #Import movie aligments from the metadata star files, which are produced by WARP at the end of processing.
        if self.doImportAlignedMovies.get():
            self.importAlignedMoviesStepId = self._insertFunctionStep('importAlignedMoviesStep',prerequisites=[self.closeSetsId])
            self.warning("Import aligned movies was set to True. Once all particles are imported, the protocol will wait for WARP to export the movie aligment star files.")

    def _stepsCheck(self):
        #We check for new steps every 5 iterations to speed up the process
        if self._stepsCheckSecsCounter < 4:
            self._stepsCheckSecsCounter += 1
            return()
        self._stepsCheckSecsCounter = 0
        #We check if the WARP importer project exists. If not, we create one.
        #If this is a continuation job, read the list of already imported micrographs.
        if not hasattr(self,'importer'):
            self._steps[0].setStatus(cons.STATUS_NEW)
            if os.path.exists(self._getExtraPath("importedMicrographs.set")):
                with open(self._getExtraPath("importedMicrographs.set"),'rb') as importedMicrographsPickleFile:
                    self._importedMicrographsSet = pickle.load(importedMicrographsPickleFile)

        newSteps = []
        newParticleSetList = self.readNewMicrographs(self.starFile.get())
        particleSetBatch = list(itertools.chain.from_iterable(newParticleSetList))
        if len(particleSetBatch):
            newSteps.append(self._insertFunctionStep('importMicrographStep',particleSetBatch,prerequisites=[]))
        if len(newSteps)>0:
            self._steps[self.closeSetsId-1].addPrerequisites(*newSteps)

        #Check if the user wants to finalize the streaming.
        close = self.streamingHasFinished()
        if close and not self._closed and self._readyToClose and not self._steps[self.closeSetsId-1].isFinished():
            self._steps[self.closeSetsId-1].setStatus(cons.STATUS_NEW)
        if hasattr(self,'importAlignedMoviesStepId'):
            if self._steps[self.importAlignedMoviesStepId-1].isRunning():
                finish = self.protocolHasFinished()
                if finish:
                    self._steps[self.importAlignedMoviesStepId-1].setStatus(cons.STATUS_FINISHED)
        self.updateSteps()

    def prepareImporterStep(self):
        '''
        This function creates the data sets that will be populated during import
        '''
        self.info("_stepsCheckSecs: {}".format(self._stepsCheckSecs))
        #self._stepsCheckSecs = max(0.1,int(self.fileTimeout.get()))
        if os.path.exists(self._getStopStreamingFilename()):
            os.remove(self._getStopStreamingFilename())

        self.importFilePath = self.starFile.get('').strip()
        if os.path.exists(self._getExtraPath("importedMicrographs.set")):
            with open(self._getExtraPath("importedMicrographs.set"),'rb') as importedMicrographsPickleFile:
                self._importedMicrographsSet = pickle.load(importedMicrographsPickleFile)

        #Creating an empty data sets
        self.outputParticles1 = SetOfParticles(filename=self._getPath("particles1.sqlite"))
        self.outputMicrographs1 = SetOfMicrographs(filename=self._getPath("micrographs1.sqlite"))
        self.outputCoordinates1 = SetOfCoordinates(filename=self._getPath("coordinates1.sqlite"))
        self.outputMovies1 = SetOfMovies(filename=self._getPath("movies1.sqlite"))
        self.outputCoordinates1.setMicrographs(self.outputMicrographs1)
        self.outputCtf1 = SetOfCTF(filename=self._getPath("ctf1.sqlite"))
        #Define outputs
        self._defineOutputs(outputParticles1 = self.outputParticles1)
        self._defineOutputs(outputMicrographs1 = self.outputMicrographs1)
        self._defineOutputs(outputCoordinates1 = self.outputCoordinates1)
        self._defineOutputs(outputCtf1 = self.outputCtf1)
        self._defineOutputs(outputMovies1 = self.outputMovies1)
        self._defineCtfRelation(self.outputMicrographs1, self.outputCtf1)
        self._defineSourceRelation(self.outputMicrographs1, self.outputCoordinates1)
        self._defineSourceRelation(self.outputMicrographs1, self.outputParticles1)
        #If movie alignments are imported, we create new, separate output datasets.
        #This is necessary in case not all movie alignments can be imported.
        #In this case, these incomplete datasets will be exposed to the user.

        self.warning("Prepared to import from {}".format(self.importFilePath))
        self.importer = WARPimporter(self, self.starFile.get(),self.outputParticles1,self.outputMicrographs1,self.outputCoordinates1, self.outputMovies1, self.outputCtf1)

    def readNewMicrographs(self, filename, **kwargs):
        '''
        This function assumes that all particles from the same microgrphs are written together in the star file.
        read from WARP goodparticles star file
        filename: The goodparticles star file
        rowToParticle: this function will be used to convert the row to Object
        '''
        newParticleSetList = []
        try:
            imgMd = md.MetaData(filename)
        except:
            print("Cant't read {}, maybe drive is busy. Skipping this iteration".format(filename))
            return(newParticleSetList)

        currentMicrograph = ""

        self.info("There are {} particles in the dataframe".format(imgMd.size()))
        for i in range(self._lastImportedId+1,imgMd.size()+1):
            row = md.Row()
            row.readFromMd(imgMd,i)
            imgRow = row.clone()
            micrographName = imgRow['rlnMicrographName']
            if micrographName not in self._importedMicrographsSet and micrographName != currentMicrograph:
                newParticleSetList.append(list())
                if currentMicrograph != "":
                    self._importedMicrographsSet.append(currentMicrograph)
                    self._lastImportedId = i
                currentMicrograph = micrographName

            if len(newParticleSetList) > 0 and micrographName == currentMicrograph:
                newParticleSetList[-1].append(imgRow)

        if self.streamingHasFinished():
            self.info("Finishing. Added {} micrographs to the import list".format(len(newParticleSetList)))
            self._importedMicrographsSet.append(currentMicrograph)
            self._readyToClose = True
            return(newParticleSetList)
        else:
            self.info("Added {} micrographs to the import list".format(len(newParticleSetList)-1))
            return(newParticleSetList[:-1]) #We always skip the last micrograph because it could be incomplete

    def importMicrographStep(self, particleSet):
        '''
        This function uses a WARPimporter object to import a set of particles
        '''
        #The import file modification time is used to decide when the import is finished.
        #In case the import file is in a different server, there might be a difference in the
        #sytem time, so using the time difference is safer.
        timeDifference = time.time() - self.mtime(self.importFilePath)
        #Create a WARP importer object with the data sets (created in prepareImporterStep())
        #Tthat will be populated by the importer
        #Start the loop
        #Import new particles
        self.importer.importParticles(particleSet)
        with open(self._getExtraPath("importedMicrographs.set"),'wb') as importedMicrographsPickleFile:
            pickle.dump(self._importedMicrographsSet, importedMicrographsPickleFile)

        #Update the output sets
        self._updateOutputSet("outputMicrographs1",self.outputMicrographs1,self.outputMicrographs1.STREAM_OPEN)
        self._updateOutputSet("outputParticles1",self.outputParticles1,self.outputParticles1.STREAM_OPEN)
        self._updateOutputSet("outputCoordinates1",self.outputCoordinates1,self.outputCoordinates1.STREAM_OPEN)
        self._updateOutputSet("outputCtf1",self.outputCtf1,self.outputCtf1.STREAM_OPEN)
        self._updateOutputSet("outputMovies1",self.outputMovies1,self.outputMovies1.STREAM_OPEN)

        #Update the summary info for the user
        summary = "Import from {} file:\n".format(self.importFilePath)

        if self.hasAttribute('outputParticles1'):
            particles = self.outputParticles1
            summary += ' Particles: *%d* ' % particles.getSize()
            summary += ('(ctf=%s, alignment=%s, phaseFlip=%s)\n'
                        % (particles.hasCTF(), particles.getAlignment(),
                           particles.isPhaseFlipped()))

        if self.hasAttribute('outputCoordinates'):
            summary += '   Coordinates: *%d* \n' % (self.outputCoordinates1.getSize())

        if self.hasAttribute('outputMicrographs1'):
            summary += '   Micrographs: *%d* \n' % (self.outputMicrographs1.getSize())

        self.summaryVar.set(summary)

    def closeSetsStep(self):
        #Before we finish this step, we update and cloaseall the data sets
        self.warning("Closing set of " + str(self.outputMicrographs1.getSize()) + "micrographs")
        self._updateOutputSet("outputMicrographs1",self.outputMicrographs1,self.outputMicrographs1.STREAM_CLOSED)
        self.warning("Closing set of " + str(self.outputParticles1.getSize()) + "particles")
        self._updateOutputSet("outputParticles1",self.outputParticles1,self.outputParticles1.STREAM_CLOSED)
        self.warning("Closing set of " + str(self.outputCoordinates1.getSize()) + "coordinates")
        self._updateOutputSet("outputCoordinates1",self.outputCoordinates1,self.outputCoordinates1.STREAM_CLOSED)
        self.warning("Closing set of " + str(self.outputCtf1.getSize()) + "CTFS")
        self._updateOutputSet("outputCtf1",self.outputCtf1,self.outputCtf1.STREAM_CLOSED)
        self.warning("Closing set of " + str(self.outputMovies1.getSize()) + "movies")
        self._updateOutputSet("outputMovies1",self.outputMovies1,self.outputMovies1.STREAM_CLOSED)
        self._closed = True

    def importAlignedMoviesStep(self):
        '''
        This function assumes that all particles from the same microgrphs are written together in the star file.
        read from WARP goodparticles star file
        filename: The goodparticles star file
        rowToParticle: this function will be used to convert the row to Object
        '''
        self.importFilePath = self.starFile.get('').strip()

        self.outputAlignedMovies1 = SetOfMovies(filename=self._getPath("alignedMovies1.sqlite"))
        self.outputParticles2 = SetOfParticles(filename=self._getPath("particles2.sqlite"))
        self.outputMicrographs2 = SetOfMicrographs(filename=self._getPath("micrographs2.sqlite"))
        self.outputCoordinates2 = SetOfCoordinates(filename=self._getPath("coordinates2.sqlite"))
        self.outputCoordinates2.setMicrographs(self.outputMicrographs2)
        self.outputCtf2 = SetOfCTF(filename=self._getPath("ctf2.sqlite"))
        self._defineOutputs(outputAlignedMovies1 = self.outputAlignedMovies1)

        particleSet = list()

        imgMd = md.MetaData(self.importFilePath)

        for i in range(1,imgMd.size()+1):
            row = md.Row()
            row.readFromMd(imgMd,i)
            imgRow = row.clone()
            particleSet.append(imgRow)

        #Create the importer object with the data sets that will be populated
        importer = WARPimporter(self, self.importFilePath, self.outputParticles2,self.outputMicrographs2,self.outputCoordinates2,self.outputAlignedMovies1,self.outputCtf2,importAlignments=True)
        #Save the time when we start waiting for the movie alingments to be available
        startTime = time.time()
        self.warning("Importing aligned movies...")
        #WARP exports movie alignment star files in the 'motion' subdirectory
        metadataPath = os.path.join(os.path.dirname(self.importFilePath), "motion")
        self.warning("Looking for metadata in {}".format(metadataPath))
        #We wait for the 'motion' subdirectory to be available
        finish = False
        while not os.path.isdir(metadataPath) and not finish:
            time.sleep(self.fileTimeout.get())
            finish = self.protocolHasFinished()
            if time.time()-startTime > self.movieTimeout.get() == 0:
                self.warning("Movie alignments not found after {}. Stop waiting. Finishing protocol".format(str(self.movieTimeout.get())))
                finish = True

        #We start impoting. It can take a while until WARP exports all movie alignments, so we do the import in a loop.
        while not finish:
            finish = self.protocolHasFinished()
            #import new particles
            importer.importParticles(particleSet)
            #Update the data sets
            if self.outputAlignedMovies1 is not None:
                self._updateOutputSet("outputAlignedMovies1",self.outputAlignedMovies1,self.outputAlignedMovies1.STREAM_OPEN)
                self._updateOutputSet("outputMicrographs2",self.outputMicrographs2,self.outputMicrographs2.STREAM_OPEN)
                self._updateOutputSet("outputParticles2",self.outputParticles2,self.outputParticles2.STREAM_OPEN)
                self._updateOutputSet("outputCtf2",self.outputCtf2,self.outputCtf2.STREAM_OPEN)
                self._updateOutputSet("outputCoordinates2",self.outputCoordinates2,self.outputCoordinates2.STREAM_OPEN)

            #If all aligned movies are available, break the loop. Else wait.
            if self.outputAlignedMovies1 is not None:
                if self.outputAlignedMovies1.getSize() == self.outputMicrographs1.getSize() or self.movieTimeout.get() == 0:
                    finish = True
                else:
                    self.warning("Not all particles have available aligned movies, waiting...")
                    time.sleep(self.fileTimeout.get())
            else:
                self.warning("Not all particles have available aligned movies, waiting...")
                time.sleep(self.fileTimeout.get())

            #If timeout is reached before all aligned movies could be imported, break the loop.
            #In this case, output incomplete data sets with available aligned movies
            if time.time()-startTime > self.movieTimeout.get():
                finish=True
                print("Movie alignments not found after {}. Stop waiting. Finishing protocol".format(str(self.movieTimeout.get())))
                #If some, but not all movies were imported, make the icomplete datasets available
                if self.outputMicrographs1.getSize() > self.outputAlignedMovies1.getSize() > 0:
                    self._defineOutputs(outputMicrographs2 = self.outputMicrographs2)
                    self._defineOutputs(outputParticles2 = self.outputParticles2)
                    self._defineOutputs(outputCoordinates2 = self.outputCoordinates2)
                    self._defineOutputs(outputCtf2 = self.outputCtf2)
                    self._defineCtfRelation(self.outputMicrographs2, self.outputCtf2)
                    self._defineSourceRelation(self.outputMicrographs2, self.outputCoordinates2)
                    self._defineSourceRelation(self.outputMicrographs2, self.outputParticles2)

        #Before we finish this step, we update and cloaseall the data sets
        self._updateOutputSet("outputAlignedMovies1",self.outputAlignedMovies1,self.outputAlignedMovies1.STREAM_CLOSED)
        self._updateOutputSet("outputMicrographs2",self.outputMicrographs2,self.outputMicrographs2.STREAM_CLOSED)
        self._updateOutputSet("outputParticles2",self.outputParticles2,self.outputParticles2.STREAM_CLOSED)
        self._updateOutputSet("outputCtf2",self.outputCtf2,self.outputCtf2.STREAM_CLOSED)
        self._updateOutputSet("outputCoordinates2",self.outputCoordinates2,self.outputCoordinates2.STREAM_CLOSED)

    # --------------------------- INFO functions -----------------------------------
    #Get the modification time of the input star file. Sometimes, this can fail if the file is
    #being updated by WARP, so we catch that exception and try again later
    def mtime(self,f):
        try:
            mtime = os.path.getmtime(f)
        except Exception as e:
            mtime = time.time()
            print(e)
            print("Star file seems to be busy. Waiting...")
        return(mtime)

    def _summary(self):
        return [self.summaryVar.get('')]

    def _methods(self):
        methods = ["Methods are not implemented"]
        return methods

    # --------------- Streaming special functions -----------------------
    def _getStopStreamingFilename(self):
        return self._getExtraPath("STOP_STREAMING.TXT")

    def _getFinishProtocolFilename(self):
        return self._getExtraPath("FINISH.TXT")

    def getActions(self):
        """ This method will allow that the 'Stop import' action to appears
        in the GUI when the user right-click in the protocol import box.
        It will allow a user to manually stop the streaming.
        """
        # Only allow to stop if running and in streaming mode
        actions = []
        if self.isRunning():
            actions.append(('STOP STREAMING', self.stopImport))
            actions.append(('FINISH NOW', self.finishProtocol))
        """
            self.info("is waiting")
            self.info(self._steps[self.closeSetsId-1].isWaiting())
            if self._steps[self.closeSetsId-1].isWaiting():

        elif self.isRunning() and hasattr(self,"importAlignedMoviesStepId"):
            if self._steps[self.importAlignedMoviesStepId-1].getStatus() == cons.STATUS_RUNNING:
        """
        return(actions)

    def stopImport(self):
        """ Since the actual protocol that is running is in a different
        we will use a simple mechanism to place an special file to stop
        process that the one that this method will be invoked from the GUI,
        the streaming.
        """
        # Just place an special file into the run folder
        self.info("Finalizing the protocol. Waiting for the last iteration to finish")
        f = open(self._getStopStreamingFilename(), 'w')
        f.close()

    def finishProtocol(self):
        f = open(self._getFinishProtocolFilename(), 'w')
        f.close()

    def streamingHasFinished(self):
        self.info("Streaming has finished: {}".format(os.path.exists(self._getStopStreamingFilename())))
        return(os.path.exists(self._getStopStreamingFilename()))

    def protocolHasFinished(self):
        self.info("Stop waiting for movie alignment files: {}".format(os.path.exists(self._getFinishProtocolFilename())))
        return(os.path.exists(self._getFinishProtocolFilename()))
