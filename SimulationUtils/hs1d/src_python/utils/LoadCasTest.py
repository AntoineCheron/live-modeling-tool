# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:31:20 2017

@author: Quentin Courtois
"""
import os
import InputsFile as IF
import OutputsFile as OF



class LoadCasTest(object):
    """
        Class reading specific Watershed test cases edited with matlab and
        saved with specific layout.

        4 files are read using this class for inputs:
            - geologic.input : f, k and d
            - hydrologic.input : t, recharge
            - key_spatialized_param.param : x_interp, Smax_interp, w_interp
            - morphologic.input : x, w, i, z_true, z_mod



        @param :
            - flag is used to select the wanted test case :
                - 0 : convergent simple watershed with synthetic recharge
                - 1 : straight simple watershed with real recharge
                - 2 : complex watershed with real recharge
                - 3 : complex watershed with synthetic recharge
                - 4 : custom watershed (needs to specify path)

            - custom_path : contains the name of the folder where is your custom
                            example
                    WARNING ! must contains only the folder's name

    """

    def __init__(self, flag=0, custom_path="", name_custom="", output=0, interp=0):

        #Setting the folder's path for the example
        test_path, self.name = self.set_path(flag, custom_path, name_custom)

        #Listing Input file's name
        self.set_inputs(test_path)

        #Listing Output file's name
        if output != 0:
            self.set_outputs(test_path)

        #Load inputs
        self.input = IF.InputsFile(self.inputs_list, interp)

        #Load outputs
        if output != 0:
            self.output = OF.OutputsFile(self.outputs_list)

    def set_path(self, flag, custom_path, name_custom):
        """
            Set the example's folder name based on user's choices
        """
        src_path = os.getcwd()
        if "\\" in src_path:
            src_path = src_path.replace("\\", "/")

        test_path = src_path[:-10] + 'test_case/matlab/'

        if flag == 0:
            test_path = test_path + '2017_1_31_8_48_45_X_conv_Y_7_slopcst_Synthetic2/'
            name = "Simple_Conv_Synth_0"
        elif flag == 1:
            test_path = test_path + '2017_1_31_9_15_18_X_straight_Y_1_slopcst_Real2/'
            name = "Simple_Straight_Real_1"
        elif flag == 2:
            test_path = test_path +  '2017_1_31_9_31_49_X_168655_Y_6784671_slopvar_Real2/'
            name = "Complex_Real_2"
        elif flag == 3:
            test_path = test_path +  '2017_1_31_9_41_29_X_362636_Y_6828317_b_slopcst_Synthetic2/'
            name = "Complex_Synth_3"
        elif flag == 4:
            if custom_path[-1] != '/':
                custom_path = custom_path + '/'
            test_path = test_path + custom_path
            name = name_custom

        return test_path, name

    def set_inputs(self, test_path):
        """
            Build a list of needed inputs files
        """
        geol_file = test_path + 'geologic.input'
        hydro_file = test_path + 'hydrologic.input'
        morpho_file = test_path + 'morphologic.input'
        spatial_param = test_path + 'key_spatialized_param.param'
        self.inputs_list = [geol_file, hydro_file, morpho_file, spatial_param]

    def set_outputs(self, test_path):
        """
            Build a list of needed outputs files
        """
        seepage_file = test_path + 'seepage_flux.code'
        storage_file = test_path + 'storage.code'
        flux_file = test_path + 'subsurface_flux.code'
        river_runoff_file = test_path + 'river_flow_runoff.output'
        relative_storage_file = test_path + 'relative_storage.state'
        self.outputs_list = [seepage_file, storage_file, flux_file, river_runoff_file, relative_storage_file]
