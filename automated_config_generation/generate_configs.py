import os
import unittest

import process_wrapper
from config import PolyploidParams


class Generate_Config_Files(unittest.TestCase):



    def test_time_series_making_configs(self):


        #sim_subfolder="sim40_1p0" #folder to make, to put put the shell scrips & qsub output
        sim_subfolder = "sim53_offdiags"
        #sim_subfolder = "sim52_diagonals"

        me_at_remote_URL='tdunn@mesx.sdsu.edu'
        template_xml_file="mesx-template.xml"
        template_sh_file="qsub-template.sh"
        local_out_dir="/home/tamsen/Data/SpecKS_input"
        commands_folder_on_mesx="/usr/scratch2/userdata2/tdunn/cluster_commands"
        output_folder_on_mesx="/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        script_destination_folder_on_mesx = os.path.join(commands_folder_on_mesx, sim_subfolder)
        specks_output_path_on_mesx = os.path.join(output_folder_on_mesx, sim_subfolder)
        origin_folder=os.path.join(local_out_dir, sim_subfolder)

        if not os.path.exists(local_out_dir):
            os.makedirs(local_out_dir)
        if not os.path.exists(origin_folder):
            os.makedirs(origin_folder)

        decimals_needed=3
        formatter = "{:0" + str(decimals_needed) + "d}"

        #for fig 1 - 4
        spec_times= [80,70, 60, 50, 40, 30, 20,10,5,0]
        #spec_times= [40]
        wgd_offsets=[0,20]
        full_sim_time  = str(100)
        #all distributions

        #for fig 5.
        #spec_times= [1,2,3,4,5,7,10,15]
        #wgd_offsets=[0,1]
        #full_sim_time=str(50)
        #dist N=5

        imp_distribution = "impulse,1,1"
        ln_distribution="lognorm,0.5,5.27"
        expon_0p1="expon,0,0.1"
        expon_1p0 = "expon,0,1"
        expon_5p0 = "expon,0,5"
        expon_10p0 = "expon,0,10"
        expon_20p0 = "expon,0,20"
        expon_100p0 = "expon,0,100"

        #coalescent_distribution_1=expon_5p0

        poly_params_by_name={}
        out_folder_by_name={}
        cluster_cmds=[]

        for i in range(0,len(spec_times)):

            job_str=str(i+1)
            for wgd_offset in wgd_offsets:

                if wgd_offset<=0:
                    job_type = "Allo"
                    #job_type = "Auto"
                else:
                    #job_type = "Allo"
                    job_type = "Auto"

                #job_type = "Allo"

                spec_time = spec_times[i]
                wgd_time = spec_time-wgd_offset

                if wgd_time < 0:
                    continue

                job_name = job_type + job_str + "_S" + formatter.format(spec_time) + "W" + formatter.format(wgd_time)
                job_params = PolyploidParams(spec_time, wgd_time,job_name)
                poly_params_by_name[job_name] = job_params
                out_folder_by_name[job_name] = os.path.join(specks_output_path_on_mesx, "specks_" + job_name)

        for poly_name,poly_params in poly_params_by_name.items():

            #coalescent_distribution= expon_1p0

            if "Allo" in poly_name: #wgd_offset < 0:, small Ne
                coalescent_distribution=  expon_5p0
            else:
                coalescent_distribution=  expon_0p1

            new_xml_file_name = poly_name + ".xml"
            xml_replacements=[("POLYPLOID_SECTION",poly_params.to_xml()),
                              ("OUTPUT_ROOT", out_folder_by_name[poly_name]),
                              ("DIV_DIST", coalescent_distribution),
                              ("FULL_SIM_TIME", full_sim_time)
                              ]
            new_file_created=write_config_file(
                template_xml_file, origin_folder,new_xml_file_name,xml_replacements)
            self.assertTrue(os.path.exists(new_file_created))
         #   / usr / scratch2 / userdata2 / tdunn / SpecKS_Output / sim9_gbd_in_auto / specks
            config_path=os.path.join(script_destination_folder_on_mesx,new_xml_file_name)
            new_script_file_name = poly_name + ".sh"
            sh_replacements=[("JOBNAME",poly_name),("CONFIGPATH",config_path)]
            new_sh_created=write_config_file(
                template_sh_file,origin_folder,new_script_file_name,sh_replacements)
            self.assertTrue(os.path.exists(new_sh_created))
            cluster_cmds.append(new_script_file_name)

        #copy over the folder of cluster commands
        cmd2=['scp','-r',origin_folder,me_at_remote_URL+':' +commands_folder_on_mesx]
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_out_dir)


        #https://stackoverflow.com/questions/26278167/submitting-a-job-to-qsub-remotely
        #https://unix.stackexchange.com/questions/8612/programmatically-creating-a-remote-directory-using-ssh
        #make the output directory
        cmd3 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;','mkdir ' + specks_output_path_on_mesx]
        print(" ".join(cmd3))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd3, local_out_dir)

        cmds_to_qsub_via_ssh=[]
        for cmd in cluster_cmds:
            cmd4 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + script_destination_folder_on_mesx + ';', 'qsub',cmd]
            cmds_to_qsub_via_ssh.append("qsub")
            cmds_to_qsub_via_ssh.append(cmd + ";")

        #submit all the qsub commands
        full_cmd_with_all_qsubs=['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + script_destination_folder_on_mesx + ';'] + cmds_to_qsub_via_ssh
        print(" ".join(full_cmd_with_all_qsubs))
        #out_string, error_string = process_wrapper.run_and_wait_with_retry(
        #    full_cmd_with_all_qsubs, local_out_dir,"Connection reset by peer",
        #    3, 10)

    # automatically set off a batch of simulation runs via qsub
    def test_making_configs(self):


        #sim_subfolder="sim40_1p0" #folder to make, to put put the shell scrips & qsub output
        #sim_subfolder = "sim50_offdiags"
        sim_subfolder = "sim60_deltat"
        me_at_remote_URL='tdunn@mesx.sdsu.edu'
        template_xml_file="mesx-template.xml"
        template_sh_file="qsub-template.sh"
        local_out_dir="/home/tamsen/Data/SpecKS_input"
        commands_folder_on_mesx="/usr/scratch2/userdata2/tdunn/cluster_commands"
        output_folder_on_mesx="/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        script_destination_folder_on_mesx = os.path.join(commands_folder_on_mesx, sim_subfolder)
        specks_output_path_on_mesx = os.path.join(output_folder_on_mesx, sim_subfolder)
        origin_folder=os.path.join(local_out_dir, sim_subfolder)

        if not os.path.exists(local_out_dir):
            os.makedirs(local_out_dir)
        if not os.path.exists(origin_folder):
            os.makedirs(origin_folder)

        decimals_needed=3
        formatter = "{:0" + str(decimals_needed) + "d}"

        #for fig 1 - 4
        #spec_times= [80,70, 60, 50, 40, 30, 20,10]
        spec_times= [40]
        wgd_offsets=[0,1,5,10,20,30,40]
        full_sim_time  = str(100)
        #all distributions

        #for fig 5.
        #spec_times= [1,2,3,4,5,7,10,15]
        #wgd_offsets=[0,1]
        #full_sim_time=str(50)
        #dist N=5

        imp_distribution = "impulse,1,1"
        ln_distribution="lognorm,0.5,5.27"
        expon_0p1="expon,0,0.1"
        expon_1p0 = "expon,0,1"
        expon_5p0 = "expon,0,5"
        expon_10p0 = "expon,0,10"
        expon_20p0 = "expon,0,20"
        expon_100p0 = "expon,0,100"

        #coalescent_distribution_1=expon_5p0

        poly_params_by_name={}
        out_folder_by_name={}
        cluster_cmds=[]

        for i in range(0,len(spec_times)):

            job_str=str(i+1)
            for wgd_offset in wgd_offsets:

                #if wgd_offset==0:
                #    job_type = "Auto"
                #else:
                #    job_type = "Allo"

                job_type = "Allo"

                spec_time = spec_times[i]
                wgd_time = spec_time-wgd_offset

                #if wgd_time < 5:
                #    continue

                job_name = job_type + job_str + "_S" + formatter.format(spec_time) + "W" + formatter.format(wgd_time)
                job_params = PolyploidParams(spec_time, wgd_time,job_name)
                poly_params_by_name[job_name] = job_params
                out_folder_by_name[job_name] = os.path.join(specks_output_path_on_mesx, "specks_" + job_name)

        for poly_name,poly_params in poly_params_by_name.items():

            coalescent_distribution= expon_1p0
            new_xml_file_name = poly_name + ".xml"
            xml_replacements=[("POLYPLOID_SECTION",poly_params.to_xml()),
                              ("OUTPUT_ROOT", out_folder_by_name[poly_name]),
                              ("DIV_DIST", coalescent_distribution),
                              ("FULL_SIM_TIME", full_sim_time)
                              ]
            new_file_created=write_config_file(
                template_xml_file, origin_folder,new_xml_file_name,xml_replacements)
            self.assertTrue(os.path.exists(new_file_created))
         #   / usr / scratch2 / userdata2 / tdunn / SpecKS_Output / sim9_gbd_in_auto / specks
            config_path=os.path.join(script_destination_folder_on_mesx,new_xml_file_name)
            new_script_file_name = poly_name + ".sh"
            sh_replacements=[("JOBNAME",poly_name),("CONFIGPATH",config_path)]
            new_sh_created=write_config_file(
                template_sh_file,origin_folder,new_script_file_name,sh_replacements)
            self.assertTrue(os.path.exists(new_sh_created))
            cluster_cmds.append(new_script_file_name)

        #copy over the folder of cluster commands
        cmd2=['scp','-r',origin_folder,me_at_remote_URL+':' +commands_folder_on_mesx]
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_out_dir)


        #https://stackoverflow.com/questions/26278167/submitting-a-job-to-qsub-remotely
        #https://unix.stackexchange.com/questions/8612/programmatically-creating-a-remote-directory-using-ssh
        #make the output directory
        cmd3 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;','mkdir ' + specks_output_path_on_mesx]
        print(" ".join(cmd3))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd3, local_out_dir)

        cmds_to_qsub_via_ssh=[]
        for cmd in cluster_cmds:
            cmd4 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + script_destination_folder_on_mesx + ';', 'qsub',cmd]
            cmds_to_qsub_via_ssh.append("qsub")
            cmds_to_qsub_via_ssh.append(cmd + ";")

        #submit all the qsub commands
        full_cmd_with_all_qsubs=['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + script_destination_folder_on_mesx + ';'] + cmds_to_qsub_via_ssh
        print(" ".join(full_cmd_with_all_qsubs))
        #out_string, error_string = process_wrapper.run_and_wait_with_retry(
        #    full_cmd_with_all_qsubs, local_out_dir,"Connection reset by peer",
        #    3, 10)


def write_config_file(template_file, out_dir, new_file_name, replacements):

    lines_to_write = []
    new_config_file = os.path.join(out_dir, new_file_name)

    with open(template_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            for r_tuple in replacements:

                out_with_the_old=r_tuple[0]
                in_with_the_new=r_tuple[1]

                if out_with_the_old in line:
                    new_line = line.replace(out_with_the_old,in_with_the_new)

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_config_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_config_file



if __name__ == '__main__':
    unittest.main()
