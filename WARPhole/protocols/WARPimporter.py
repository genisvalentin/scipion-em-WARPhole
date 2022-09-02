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

import os
from collections import OrderedDict
from emtable import Table
import time
import copy

from pyworkflow.object import Float
from pwem.constants import ALIGN_PROJ, ALIGN_2D, ALIGN_NONE
from pwem.objects import Micrograph, MovieAlignment, Movie, Particle, Coordinate
import pwem.emlib.metadata as md
import pyworkflow.utils as pwutils
from relion.convert.convert_deprecated import setupCTF

from relion.convert.convert_utils import relionToLocation, locationToRelion
from relion.convert.convert_deprecated import rowToParticle, rowToCoordinate, rowToCtfModel
from relion.convert.convert31 import OpticsGroups


'''
Need to fix the following errors
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMovies1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMovies1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMovies1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMovies1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMovies1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputMovies1 has no sampling rate!!!
FATAL ERROR: Object 1641.outputParticles1 has no sampling rate!!!
'''

'''
Something to fix!
FATAL ERROR: Object 5750.outputParticles1 has no sampling rate!!!
FATAL ERROR: Object 5750.outputMicrographs1 has no sampling rate!!!
FATAL ERROR: Object 5750.outputMovies1 has no sampling rate!!!

Something to fix!
Check for new files after the correct fileTimeout
I think I do not need the 0 = no streaming option anymore. I have now a different
mechanism to stop streaming.
'''

class WARPimporter:
    """ Helper class to import WARP-generated particles in streaming mode """
    def __init__(self, protocol, starFile, partSet, micSet=None, coordSet=None, movieSet=None, ctfSet=None, importAlignments=False):
        self.protocol = protocol
        self._starFile = starFile
        self.copyOrLink = self.protocol.copyBinaries.get()
        self.version30 = False
        self._imgDict = {}
        self._micDict = {}
        self._movieDict = {}
        self._coordDict = {}
        self.partSet = partSet
        self.micSet = micSet
        self.movieSet = movieSet
        self.coordSet = coordSet
        self.ctfSet = ctfSet
        self._importedMovies = set()
        self._importedParticles = set()
        self._importedMicrographs = set()
        self._importAlignments = importAlignments
        self._importedCoords = set()
        self.acqRow = None
        self.range = None
        self.dim = None
        self._initSets()

    def _initSets(self):
        '''This function prepares the particle, coordinate, movie and micrographs sets'''
        if self.acqRow is None:
            _, self.modelRow, self.acqRow = self._findImagesPath('rlnImageName')

        if self.micSet is not None:
            self.micSet.setObjComment('Averaged micrographs imported from Relion star file:\n%s' % self._starFile)
            self.micSet.enableAppend()
            if not hasattr(self,'acquisitionDict'):
                self.loadAcquisitionInfo(self.micSet)
            if self.micSet.getSize() > 0:
                self.micSet.loadAllProperties()
            else:
                self.loadAcquisitionInfo(self.micSet)
                self.micSet.setSamplingRate(self.acquisitionDict['micrographSamplingRate'])

        if self.partSet is not None:
            self.partSet.setObjComment('Particles imported from Relion star file:\n%s' % self._starFile)
            self.partSet.enableAppend()
            if self.partSet.getSize() > 0:
                self.partSet.loadAllProperties()
            else:
                self.loadAcquisitionInfo(self.partSet)
                self.partSet.setSamplingRate(self.acquisitionDict['particleSamplingRate'])

        if self.movieSet is not None:
            self.movieSet.setObjComment('Movies imported from Relion star file:\n%s' % self._starFile)
            self.movieSet.enableAppend()
            if self.movieSet.getSize() > 0:
                self.movieSet.loadAllProperties()
            else:
                self.movieSet.setSamplingRate(self.acquisitionDict['movieSamplingRate'])
                self.loadAcquisitionInfo(self.movieSet)

        if self.coordSet is not None:
            self.coordSet.setObjComment('Particle coordinates imported from Relion star file:\n%s' % self._starFile)
            self.coordSet.enableAppend()
            if self.coordSet.getSize() > 0:
                self.coordSet.loadAllProperties()

        if self.ctfSet is not None:
            self.ctfSet.setObjComment('CTF models imported from Relion star file:\n%s' % self._starFile)
            self.ctfSet.enableAppend()
            if self.ctfSet.getSize() > 0:
                self.ctfSet.loadAllProperties()

    def importParticles(self,particleSet):
        '''Main method of this class. Needs to be called to import particles'''
        self._initSets()
        newFiles = self.readSetOfNewParticles(
            particleSet,
            preprocessImageRow=self._preprocessImageRow30,
            postprocessImageRow=self._postprocessImageRow30,
            readAcquisition=False)
        if self.coordSet is not None and self.partSet is not None:
            #self.coordSet.setBoxSize(self.partSet.getDimensions()[0])
            self.coordSet.setBoxSize(self.partSet.getDimensions()[0]*round(self.acquisitionDict['particleSamplingRate']/self.acquisitionDict['micrographSamplingRate']))
        self._importedParticles = newFiles
        self.protocol.info("Added {} new particles".format(str(len(newFiles))))

    def _findImagesPath(self, label, warnings=True):
        '''This function validates the input path for the binaries and gets the acquisition settings from the first row'''
        # read the first table
        try:
            table = Table(fileName=self._starFile, tableName='particles')
        except:
            table = Table(fileName=self._starFile)

        acqRow = row = table[0]
        if row is None:
            raise Exception("Cannot import from empty metadata: %s"
                            % self._starFile)

        if not row.get('rlnOpticsGroup', False):
            self.version30 = True
            self.protocol.warning("Import from Relion version < 3.1 ...")
        else:
            # acqRow = OpticsGroups.fromStar(self._starFile)
            # read particles table
            table = Table(fileName=self._starFile, tableName='particles')
            row = table[0]

        self.protocol.info("Reading row")
        self.protocol.info(row)

        if not row.get(label, False):
            raise Exception("Label *%s* is missing in metadata: %s"
                            % (label, self._starFile))

        index, fn = relionToLocation(row.get(label))
        self._imgPath = pwutils.findRootFrom(self._starFile, fn)
        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s"
                                  % self._starFile)
        # Check if the MetaData contains either MDL_MICROGRAPH_ID
        # or MDL_MICROGRAPH, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = (row.get('rlnMicrographName', False) or
                             row.get('rlnMicrographId', False))

        print("acqRow",acqRow)

        index, imgPath = relionToLocation(table[0].get('rlnMicrographName'))
        print("Reading movie dimensions from {}".format(os.path.join(self._imgPath, imgPath)))
        m = Movie(os.path.join(self._imgPath, imgPath))
        self.dim = m.getDim()
        print("Dim: (%s)" % ", ".join(map(str, self.dim)))
        self.range = [1, self.dim[2], 1]
        self.movieSet.setDim(self.dim)
        self.movieSet.setFramesRange(self.range)

        return row, None, acqRow

        #This function imports the movies and the micrographs
    def _preprocessImageRow30(self, img, imgRow):
        self.preprocess_success = False

        #Create a link or copy the particle binary files (*.mrcs). If the file already exists, does nothing
        self.copyOrLinkBinary(imgRow, 'rlnImageName', self._imgPath, self.protocol._getExtraPath(), copyFiles=self.protocol.copyBinaries.get())
        setupCTF(imgRow, self.acquisitionDict['micrographSamplingRate'])

        movieId = imgRow.get('rlnMicrographId', None)
        movieName = imgRow.get('rlnMicrographName', None)

        #Import the associated movie (*.mrcs)
        if self.movieSet is not None and imgRow.get('rlnMicrographName', None) is not None:
            # Check which is the key to identify micrographs (id or name)
            if movieId is not None:
                movieKey = movieId
            else:
                movieKey = movieName

            # First time I found this micrograph (either by id or name)
            if movieKey not in self._importedMovies:
                newName = self.copyOrLinkBinary(imgRow, 'rlnMicrographName', self._imgPath, self.protocol._getExtraPath(), copyFiles=self.protocol.copyBinaries.get())
                movieName = imgRow.get('rlnMicrographName', None)
                movie = Movie()
                movie.setObjId(movieId)
                movie.setFileName(newName)
                movie.setMicName(newName)
                if self.range:
                    movie.setFramesRange(self.range)

                if self._importAlignments:
                    alignment = self.getMicrographAlignment(movie)
                    if alignment:
                        movie.setAlignment(alignment)
                        self.movieSet.append(movie)
                        # Update dict with new movie
                        self._movieDict[movieKey] = movie
                        self._importedMovies.add(movieKey)
                    else:
                        self.preprocess_success = False
                        return()
                else:
                    self.movieSet.append(movie)
                    # Update dict with new movie
                    self._movieDict[movieKey] = movie
                    self._importedMovies.add(movieKey)


        if self.micSet is not None:
            imgRow['rlnMicrographName'] = self.fixMicName(imgRow['rlnMicrographName'])
            micName = imgRow.get('rlnMicrographName', None)
            micId = imgRow.get('rlnMicrographId', None)
            # Check which is the key to identify micrographs (id or name)
            if micId is not None:
                micKey = micId
            else:
                micKey = micName

            # First time I found this micrograph (either by id or name)
            if micKey not in self._importedMicrographs:
                self.copyOrLinkBinary(imgRow, 'rlnMicrographName', self._imgPath, self.protocol._getExtraPath(), copyFiles=self.protocol.copyBinaries.get())
                micName = imgRow.get('rlnMicrographName', None)
                self.protocol.info("Importing from micrograph {}".format(micName))
                mic = Micrograph()
                mic.setObjId(micId)
                if micName is None:
                    micName = self.protocol._getExtraPath('fake_micrograph%6d' % micId)
                mic.setFileName(self.protocol._getExtraPath(micName))
                mic.setMicName(movieName)
                ctf = rowToCtfModel(imgRow)
                ctf.setMicrograph(mic)
                self.ctfSet.append(ctf)
                mic.setCTF(rowToCtfModel(imgRow))
                self.micSet.append(mic)
                self._micDict[os.path.basename(movieName)] = mic
                self._importedMicrographs.add(micKey)
            else:
                mic = self._micDict.get(os.path.basename(movieName))
                movieName = mic.getMicName()

            # Update the row to set a MDL_MICROGRAPH_ID
            imgRow['rlnMicrographId'] = int(mic.getObjId())
            imgRow['rlnMicrographName'] = movieName
            img.setCTF(rowToCtfModel(imgRow))
            self.preprocess_success = True

    #This function imports the particle coordinates
    def _postprocessImageRow30(self, img, imgRow):
        if self._micIdOrName is not None and self.coordSet is not None:
            micId = imgRow.get('rlnMicrographId', None)
            micName = imgRow.get('rlnMicrographName', None)
            partName = imgRow.get('rlnImageName',None)
            if img.hasCoordinate():
                coord = img.getCoordinate()
                #Adjust the coordinate in pixels in case particles are binned. Rounding because binnig is always an integer.
                coord.setMicId(micId)
                coord.setMicName(micName)
                if partName not in self._importedCoords and self.preprocess_success:
                    self._importedCoords.add(partName)
                    scaled_coord = Coordinate()
                    scaled_coord.copyInfo(coord)
                    scaled_coord.setMicId(micId)
                    scaled_coord.setMicName(micName)
                    scaled_coord.scale(round(self.acquisitionDict['particleSamplingRate']/self.acquisitionDict['micrographSamplingRate']))
                    self._coordDict[partName] = scaled_coord
                    self.coordSet.append(scaled_coord)

    #Return a dictionary with acquisition values and the sampling rate information.
    #This informatoin is taken from the first particle of th star file.
    def loadAcquisitionInfo(self,micSet):
        acquisitionDict = {}
        print("Getting acquisition info")

        try:
            #_, modelRow, acqRow = self._findImagesPath('rlnImageName')
            acqRow = self.acqRow
            #print("acqRow",acqRow)
            acquisition = micSet.getAcquisition()
            if acqRow.get('rlnVoltage', False):
                acquisitionDict['voltage'] = acqRow.rlnVoltage
                acquisition.setVoltage(acquisitionDict['voltage'])

            if acqRow.get('rlnAmplitudeContrast', False):
                acquisitionDict['amplitudeContrast'] = acqRow.rlnAmplitudeContrast
                acquisition.setAmplitudeContrast(acquisitionDict['amplitudeContrast'])

            if acqRow.get('rlnSphericalAberration', False):
                acquisitionDict['sphericalAberration'] = acqRow.rlnSphericalAberration
                acquisition.setSphericalAberration(acquisitionDict['sphericalAberration'])

            if acqRow.get('rlnDetectorPixelSize', False):
                acquisitionDict['particleSamplingRate'] = acqRow.rlnDetectorPixelSize


            acquisitionDict['movieSamplingRate'] = self.protocol.moviePixelSize.get()
            acquisitionDict['micrographSamplingRate'] = self.protocol.micrographPixelSize.get()
            if acquisitionDict['micrographSamplingRate'] == -1:
                acquisitionDict['micrographSamplingRate'] = acquisitionDict['particleSamplingRate']

            acquisition.setDosePerFrame(self.protocol.dosePerFrame.get())
            acquisition.setMagnification(self.protocol.magnification.get())
            micSet.setAcquisition(acquisition)

        except Exception as ex:
            print("Error loading acquisition: ", str(ex))

        self.acquisitionDict = acquisitionDict

    #Returns the path to the micrograph metadata. We assume that the metadata star files are located in motion/*.star
    def getMicrographMetadata(self,micName):
        metadataName = pwutils.replaceBaseExt(micName, 'star')
        metadataPath = os.path.join(os.path.dirname(self._starFile), "motion", metadataName)
        print("Loading metadata file {}".format(metadataPath))
        if not os.path.isfile(metadataPath):
            print("Metadata for {} not found".format(micName))
            return(False)
        return(metadataPath)

    #Reads the micrograph metadata and returns the movie alignment
    def getMicrographAlignment(self,movie):
        movieName = movie.getMicName()
        motionStar = self.getMicrographMetadata(movieName)
        if os.path.isfile(motionStar):
            table = Table(fileName=motionStar, tableName='global_shift')
            xshifts, yshifts, frames = [], [], []
            for row in table:
                frames.append(row.rlnMicrographFrameNumber)
                xshifts.append(row.rlnMicrographShiftX)
                yshifts.append(row.rlnMicrographShiftY)
            alignment = MovieAlignment(first=frames[0], last=frames[-1], xshifts=xshifts, yshifts=yshifts)
            alignment.setRoi([0,0,0,0])
            return(alignment)
        else:
            return(False)

    #Fixes the rlnMicrographName from tiff into average/*.mrc
    def fixMicName(self,micName):
        return(os.path.join("average", pwutils.replaceBaseExt(micName, 'mrc')))

    #Reads the goodparticles star file generated by WARP and imports the new particles
    def readSetOfNewParticles(self, particleSet, **kwargs):
        """read from WARP goodparticles star file
            filename: The goodparticles star file
            rowToParticle: this function will be used to convert the row to Object
        """
        img = None
        newFiles = set()
        oldFiles = set(self._imgDict.keys())
        for imgRow in particleSet:
            imgName = imgRow['rlnImageName']
            #img = self._imgDict.get(imgName, None)
            if imgName not in oldFiles:
                img = rowToParticle(imgRow, **kwargs)
                if not self.preprocess_success:
                    continue
                self._imgDict[imgName] = img
                newFiles.add(imgName)
                self.partSet.append(img)

        if not img is None:
            self.partSet.setHasCTF(img.hasCTF())

        return(newFiles)

    #Create a symlink or copy the binary files (particles, micrographs, movies) into the Scipion project dir
    def copyOrLinkBinary(self, imgRow, label, basePath, destBasePath ,copyFiles=False):
        index, imgPath = relionToLocation(imgRow.get(label))
        baseName = os.path.join(os.path.dirname(imgPath),os.path.basename(imgPath))
        os.makedirs(os.path.join(destBasePath,os.path.dirname(imgPath)),exist_ok=True)
        newName = os.path.join(destBasePath, baseName)
        if not os.path.exists(newName):
            if copyFiles:
                pwutils.copyFile(os.path.join(basePath, imgPath), newName)
            else:
                pwutils.createLink(os.path.join(basePath, imgPath), newName)
        if not label=="rlnMicrographName":
            imgRow.set(label, locationToRelion(index, newName))
        return(newName)
