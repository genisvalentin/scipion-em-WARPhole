"""
Create an xml_dir
Launch the script to split into optic groups every x micrographs
The script generates a star files
Parse the star file
Keep a dictionary of micrographs and optics groups that are parsed from the star file
update the optic groups number (add depending on the number of iterations)
Use the star file to add the optic groups to the open particle set (create a new copy)
Re-use the code from:
https://scipion-em.github.io/docs/_modules/relion/protocols/protocol_assign_optic_groups.html#ProtRelionAssignOpticsGroup
"""

import pyworkflow.utils as pwutils
import os
from datetime import datetime
import time
import emtable
from pyworkflow.object import Integer
from relion.convert.convert31 import OpticsGroups, getPixelSizeLabel

import pyworkflow.protocol.constants as cons
from pyworkflow import VERSION_2_0
from pwem.protocols import EMProtocol
from pyworkflow.object import Set
from pyworkflow.protocol.params import BooleanParam, IntParam, PointerParam, GT, FolderParam
from xmipp3.protocols.protocol_trigger_data import XmippProtTriggerData

from .EPU_Group_AFIS import main as EPU_Group_AFIS_main

class AssignOpticsGroup(XmippProtTriggerData):
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
    _label = 'assign optic groups'
    _lastUpdateVersion = VERSION_2_0
    micDict = {'': 0}

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('inputImages', PointerParam,
                      pointerClass='SetOfImages',
                      label='Input images', important=True)

        form.addParam('XMLpath', FolderParam, label="Path to XML files", important=True)

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

        micDict = self.assignEPUGroupAFIS(self.newImages,str(self.XMLpath))
        micDict =self. shiftMicDict(micDict,max(self.micDict.values()))
        self.micDict = {**self.micDict, **micDict}
        self.addOpticsGroup(self.newImages,self.micDict)

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

    #################### Utility fucntions #####################
    def assignEPUGroupAFIS(self,partSet,XMLpath):
        subfolder = self.importXmlFiles(partSet,XMLpath)
        micDict = {}
        if subfolder:
            starFile = os.path.join(subfolder,"optics.star")
            self.info("Running script in subfolder {} and stafile {}".format(subfolder,starFile))
            self.runAFISscript(subfolder, starFile)
            micDict = self.readOpticsGroupStarFile(starFile)
        return(micDict)

    def importXmlFiles(self,partSet,XMLpath):
        self.info("Looking for XML files in {}".format(XMLpath))
        partPaths = [part.getFileName() for part in partSet]
        XMLpaths = set([pwutils.path.replaceBaseExt(p, 'xml') for p in partPaths])
        subfolder = self._getExtraPath(datetime.now().strftime("%H%M%S"))
        for path in XMLpaths:
            pwutils.path.makeFilePath(os.path.join(subfolder,path))
        counter = 0
        while True: ##Wait until all XML files are available, or give after 10 attempts
            wait = False
            counter += 1
            for p in XMLpaths:
                if not os.path.isfile(os.path.join(XMLpath,p)):
                    wait = True
            if wait and counter < 11:
                self.info("Waiting for XML files, sleeping for 60 seconds")
                time.sleep(self.delay)
            else:
                break
        self.info("Importing XML files")
        return_val = 0
        for p in XMLpaths:
            if os.path.isfile(os.path.join(XMLpath,p)):
                pwutils.path.copyFile(os.path.join(XMLpath,p), os.path.join(subfolder,p))
                return_val = subfolder
        return(return_val)

    def runAFISscript(self,XMLpath,outputStarFile):
        self.info("Running EPU_Group_AFIS script")
        EPU_Group_AFIS_main(xml_dir = XMLpath, output_fn = outputStarFile)

    def readOpticsGroupStarFile(self,starFile):
        self.info("Reading optics groups from file: {}".format(starFile))
        og = OpticsGroups.fromStar(starFile)
        micTable = emtable.Table(fileName=inputStar,tableName='micrographs')
        micDict = {row.rlnMicrographName: row.rlnOpticsGroup for row in micTable}
        return(micDict)

    def shiftMicDict(self,micDict,i):
        for key, value in micDict.items():
            micDict[key] = value + i
        return(micDict)

    def addOpticsGroup(self,partSet,micDict):
        self.info("Updating optics groups in output particle set")
        for mic in partSet:
            ogNumber = micDict.get(mic.getFileName(),1)
            if not hasattr(mic, '_rlnOpticsGroup'):
                mic._rlnOpticsGroup = Integer()
            mic._rlnOpticsGroup.set(Integer(ogNumber))
