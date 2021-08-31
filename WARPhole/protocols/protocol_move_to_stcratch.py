'''
create an empty particle stack
read open particle set and get at list of new stacks
if new stacks, iterate through the list of new stacks and insertfunctionstep to move each stack
	get the following particle stack, if there is none, go back to sleep
	check if there is enough disk space for the destination stacks
		if not enough space, stop here
	put the particles in the stack into a loading port
	spawn a new copy operation for the next stack
	when copy finished, make the new particle stack available
else
	sleep
go back to 2
'''
# **************************************************************************
# *
# * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
# *              David Maluenda (dmaluenda@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import time
from datetime import datetime

import pyworkflow.protocol.constants as cons
from pyworkflow import VERSION_2_0
from pwem.protocols import EMProtocol
from pyworkflow.object import Set
from pyworkflow.protocol.params import BooleanParam, IntParam, PointerParam, GT

class MoveToScratch(EMProtocol):
    """
    Waits until certain number of images is prepared and then
    send them to output.
    It can be done in 3 ways:
        - If *Send all particles to output?*' is _No_:
            Once the number of images is reached, a setOfImages is returned and
            the protocols finished (ending the streaming from this point).
        - If *Send all particles to output?*' is _Yes_ and:
            - If *Split particles to multiple sets?* is _Yes_:
                Multiple closed outputs will be returned as soon as
                the number of images is reached.
            - If *Split particles to multiple sets?* is _No_:
                Only one output is returned and it is growing up in batches of
                a certain number of images (completely in streaming).
    """
    _label = 'trigger data'
    _lastUpdateVersion = VERSION_2_0

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('inputImages', PointerParam,
                      pointerClass='SetOfImages',
                      label='Input images', important=True)
        form.addParam('outputSize', IntParam, default=10000,
                      label='Minimum output size',
                      help='How many particles need to be on input to '
                           'create output set.')
        form.addParam('allImages', BooleanParam, default=True,
                      label='Send all items to output?',
                      help='If NO is selected, only a closed subset of '
                           '"Output size" items will be send to output.\n'
                           'If YES is selected it will still running in streaming.')
        form.addParam('splitImages', BooleanParam, default=False,
                      label='Split items to multiple sets?',
                      condition='allImages',
                      help='If YES is selected, multiple closed outputs of '
                           '"Output size" are returned.\n'
                           'If NO is selected, only one open and growing output '
                           'is returned')
        form.addParam('delay', IntParam, default=10, label="Delay (sec)",
                      validators=[GT(3, "must be larger than 3sec.")],
                      help="Delay in seconds before checking new output")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # initializing variables
        self.finished = False
        self.images = []
        self.splitedImages = []
        self.outputCount = 0
        self.setImagesClass()
        self.setImagesType()

        # steps
        imsSteps = self._insertFunctionStep('delayStep')
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=[imsSteps], wait=True)

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def createOutputStep(self):
        pass


    def _checkNewInput(self):
        imsFile = self.inputImages.get().getFileName()
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(imsFile))

        # If the input's sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'newImages'):
            return None

        # loading the input set in a dynamic way
        inputClass = self.getImagesClass()
        self.imsSet = inputClass(filename=imsFile)
        self.imsSet.loadAllProperties()

        # loading new images to process
        if len(self.images) > 0:  # taking the non-processed yet
            self.newImages = [m.clone() for m in self.imsSet.iterItems(
                orderBy='creation',
                where='creation>"' + str(self.check) + '"')]
        else:  # first time
            self.newImages = [m.clone() for m in self.imsSet]
        self.moveImages(self.newImages)
        self.splitedImages = self.splitedImages + self.newImages
        self.images = self.images + self.newImages
        if len(self.newImages) > 0:
            for item in self.imsSet.iterItems(orderBy='creation',
                                              direction='DESC'):
                self.check = item.getObjCreation()
                break

        self.lastCheck = datetime.now()
        self.streamClosed = self.imsSet.isStreamClosed()
        self.imsSet.close()

        # filling the output if needed
        self._fillingOutput()

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        if self.streamClosed:
            self.finished = True
        elif not self.allImages.get():
            self.finished = len(self.images) >= self.outputSize
        else:
            self.finished = False

        outputStep = self._getFirstJoinStep()
        deps = []
        if self.finished:  # Unlock createOutputStep if finished all jobs
            self._fillingOutput()  # To do the last filling
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)
        else:
            delayId = self._insertFunctionStep('delayStep', prerequisites=[])
            deps.append(delayId)

        if outputStep is not None:
            outputStep.addPrerequisites(*deps)
        self.updateSteps()

    def _fillingOutput(self):
        imsSqliteFn = '%s.sqlite' % self.getImagesType('lower')
        outputName = self.getOututName()
        if len(self.images) >= self.outputSize or self.finished:
            if self.allImages:  # Streaming and semi-streaming
                if self.splitImages:  # Semi-streaming: Splitting the input
                    if len(self.splitedImages) >= self.outputSize or \
                            (self.finished and len(self.splitedImages) > 0):
                        self.outputCount += 1
                        imageSet = self._loadOutputSet(self.getImagesClass(),
                                                       '%s%d.sqlite'
                                                       % (self.getImagesType('lower'),
                                                          self.outputCount),
                                                       self.splitedImages)
                        # The splitted outputSets are always closed
                        self._updateOutputSet("%s%d" % (outputName, self.outputCount),
                                              imageSet, Set.STREAM_CLOSED)
                        self.splitedImages = []
                else:  # Full streaming case
                    if not os.path.exists(self._getPath(imsSqliteFn)):
                        imageSet = self._loadOutputSet(self.getImagesClass(),
                                                       imsSqliteFn,
                                                       self.images)
                    else:
                        # if finished no images to add, but we need to close the set
                        imagesToAdd = self.newImages if not self.finished else []
                        imageSet = self._loadOutputSet(self.getImagesClass(),
                                                       imsSqliteFn,
                                                       imagesToAdd)
                    streamMode = Set.STREAM_CLOSED if self.finished else \
                        Set.STREAM_OPEN
                    self._updateOutputSet(outputName, imageSet, streamMode)

            elif not os.path.exists(self._getPath(imsSqliteFn)):
                imageSet = self._loadOutputSet(self.getImagesClass(), imsSqliteFn,
                                               self.images)
                # The outputSet is always closed here
                self._updateOutputSet(outputName, imageSet, Set.STREAM_CLOSED)

    def _loadOutputSet(self, SetClass, baseName, newImages):
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputImages.get()
        outputSet.copyInfo(inputs)
        outputSet.copyItems(newImages)
        return outputSet

    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
        outputSet.setStreamState(state)
        if self.hasAttribute(outputName):
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protocol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            self._defineOutputs(**{outputName: outputSet})
            self._defineTransformRelation(self.inputImages, outputSet)
            self._store(outputSet)

        # Close set database to avoid locking it
        outputSet.close()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        outputStr = self.getOututName()

        if self.allImages.get() and not self.splitImages.get():
            summary.append("MODE: *full streaming*.")
            triggeredMsg = ("'%s' released, it will be growing up "
                            "as soon as the input does." % outputStr)
        elif self.splitImages.get():
            summary.append("MODE: *semi streaming (batches)*.")
            triggeredMsg = ("%d '%s' are being released with %d items, "
                            "each. A new batch will be created when ready."
                            % (self.getOutputsSize(), outputStr, self.outputSize))
        else:
            summary.append("MODE: *static output*.")
            triggeredMsg = ("'%s' released and closed. Nothing else to do."
                            % outputStr)

        if self.getOutputsSize():
            summary.append(triggeredMsg)
        else:
            inputStr = self.getImagesType() if self.inputImages.get() else '(or not ready)'
            summary.append("Not enough input %s to release an output, yet."
                           % inputStr)
            summary.append("At least, %d items are needed to trigger an output."
                           % self.outputSize.get())

        if (self.isFinished() and
                self.outputSize.get() > [o for o in
                                         self.iterOutputAttributes()][0][1].getSize()):
            summary.append("Output released because streaming finished.")

        return summary

    def _validate(self):
        pass

    def _moveImages(imgSet):
        imgSetSize = self._getImgSetSize(imtSet)
        freeScratchSpace = self._getFreeScratchSpace(imgSet)
        print("imgSetSize: {}, freeScratchSpace: {}".format(str(imgSetSize),str(freeScratchSpace)))
        while imgSetSize > freeScratchSpace:
            time.sleep(60)
            freeScratchSpace = self._getFreeScratchSpace()
            print("Not enough scratch space available. Sleeping for 60 seconds")
        for img in imgSet:
            filename = img.getFileName()
            path,name = os.path.dirname(filename),os.path.basename(filename)
            newFilename = os.path.join(scratchPath,filename)
            pyworkflow.utils.path.makeFilePath(newFilename)
            symlink = self.protocol._getExtraPath(filename)
            pyworkflow.utils.path.makeFilePath(symlink)
            if not os.path.exists(newFilename):
            	print("Moving {} to {}".format(filename,newFilename))
            	pyworkflow.utils.path.copyFile(filename, newFilename)
            if not os.path.exists(symlink):
            	print("Creating symling from {} to {}".format(symlink,newFilename))
            	pyworkflow.utils.path.createLink(symlink, newFilename)
            img.setFileName(symlink)

    def _getImgSetSize(imgSet):
        totalSize = 0
        for img in imgSet:
            totalSize += pyworkflow.utils.path.getFileSize(img.getFileName())
        return(totalSize)

    def _getFreeScratchSpace(scratchPath):
        return(shutil.disk_usage(path).free)

    # --------------------------- UTILS functions -----------------------------
    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def delayStep(self):
        time.sleep(self.delay)

    def setImagesClass(self):
        self._inputClass = self.inputImages.get().getClass()

    def getImagesClass(self):
        return self._inputClass

    def setImagesType(self):
        inputSet = self.inputImages.get()
        if inputSet:
            inputClassName = inputSet.getClassName()
            self._inputType = inputClassName.split('SetOf')[1]
        else:
            self._inputType = None

    def getImagesType(self, letters='default'):
        if not self.hasAttribute('_inputType'):
            self.setImagesType()
        typeStr = str(self._inputType)
        if letters == 'lower':
            return typeStr.lower()
        else:
            return typeStr

    def getOututName(self):
        return 'output%s' % self.getImagesType()
