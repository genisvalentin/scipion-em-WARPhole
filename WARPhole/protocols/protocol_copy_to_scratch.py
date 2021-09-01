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
import time
from datetime import datetime
import shutil

import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
from pyworkflow import VERSION_2_0
from pwem.protocols import EMProtocol
from pyworkflow.object import Set
from pyworkflow.protocol.params import BooleanParam, IntParam, PointerParam, GT, FolderParam
from xmipp3.protocols.protocol_trigger_data import XmippProtTriggerData

class CopyToScratch(XmippProtTriggerData):
    """
	Moves particle mrcs files to the scratch drive.
	The output particle sets are modified so that the particle filenames are
	symlinks pointing to the scratch drive data.
	If revert is "Yes", one can revert the symlinks to the original filenames.
	This is useful if the scratch drive data is not available anymore.
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
    _label = 'copy to scratch'
    _lastUpdateVersion = VERSION_2_0

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('inputImages', PointerParam,
                      pointerClass='SetOfImages',
                      label='Input images', important=True)

        form.addParam('revert',
		                      BooleanParam,
		                      label='Revert',
		                      default=False,
		                      help="Whether to revert the particle set from the scratch drive to the original directory")

        form.addParam('scratchPath', FolderParam, label="Scratch directory", important=True, condition='(revert == False)')

        form.addParam('outputSize', IntParam, default=10000, condition='(revert == False)',
                      label='Minimum output size',
                      help='How many particles need to be on input to '
                           'create output set.')
        form.addParam('allImages', BooleanParam, default=True, condition='(revert == False)',
                      label='Send all items to output?',
                      help='If NO is selected, only a closed subset of '
                           '"Output size" items will be send to output.\n'
                           'If YES is selected it will still running in streaming.')
        form.addParam('splitImages', BooleanParam, default=False,
                      label='Split items to multiple sets?',
                      condition='allImages and revert == False',
                      help='If YES is selected, multiple closed outputs of '
                           '"Output size" are returned.\n'
                           'If NO is selected, only one open and growing output '
                           'is returned')
        form.addParam('delay', IntParam, default=10, label="Delay (sec)", condition='(revert == False)',
                      validators=[GT(3, "must be larger than 3sec.")],
                      help="Delay in seconds before checking new output")

    # --------------------------- INSERT steps functions ----------------------

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
        if self.revert:
            self._revertImages(self.newImages)
        else:
            self._moveImages(self.newImages)
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

    def _moveImages(self,imgSet):
        scratchPath = str(self.scratchPath)
        imgSetSize = self._getImgSetSize(imgSet)
        freeScratchSpace = self._getFreeScratchSpace(scratchPath)
        print("imgSetSize: {}, freeScratchSpace: {}".format(str(imgSetSize),str(freeScratchSpace)))
        while imgSetSize > freeScratchSpace:
            time.sleep(60)
            freeScratchSpace = self._getFreeScratchSpace()
            print("Not enough scratch space available. Sleeping for 60 seconds")
        for img in imgSet:
            filename = img.getFileName()
            newFilename = os.path.join(scratchPath,filename)
            pwutils.path.makeFilePath(newFilename)
            symlink = self._getExtraPath(filename)
            pwutils.path.makeFilePath(symlink)
            if not os.path.exists(newFilename):
            	print("Copying {} to {}".format(filename,newFilename))
            	pwutils.path.copyFile(filename, newFilename)
            if not os.path.exists(symlink):
            	print("Creating symlink from {} to {}".format(symlink,newFilename))
            	pwutils.path.createLink(newFilename, symlink)
            img.setFileName(symlink)

    def _revertImages(self,imgSet):
        for img in imgSet:
            filename = img.getFileName()
            newFilename = 'Runs'.join([filename.split("Runs")[0]] + filename.split("Runs")[2:])
            img.setFileName(newFilename)

    def _getImgSetSize(self,imgSet):
        totalSize = 0
        for img in imgSet:
            totalSize += pwutils.path.getFileSize(img.getFileName())
        return(totalSize)

    def _getFreeScratchSpace(self,path):
        return(shutil.disk_usage(path).free)
