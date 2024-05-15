import os
import time
import unittest
import process_wrapper
from results_viewer.download_sim_results import get_subdirectories, get_run_folders_by_polyploid_name


class MyTestDownloaderForFits(unittest.TestCase):



    def test_download_mesx_results(self):

        polyploid_names = ['Allo_Coffea15','Allo_Poplar16','Allo14_Maize']
        batch_folders = ["sim41_coffee", "sim41_Poplar","sim41_Maize"]
        me_at_remote_URL =  'tdunn@mesx.sdsu.edu'
        local_output_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        remote_output_folder = "/usr/scratch2/userdata2/tdunn/SpecKS_Output"

        for i in range(0,len(batch_folders)):
            batch_folder=batch_folders[i]
            polyploid_name=polyploid_names[i]
            local_batch_folder = os.path.join(local_output_folder, batch_folder)
            remote_batch_folder = os.path.join(remote_output_folder, batch_folder)
            get_SpecKS_output_for_fit_test(local_batch_folder, local_output_folder, me_at_remote_URL, polyploid_name,
                                            remote_batch_folder)

def get_SpecKS_output_for_fit_test(local_batch_folder, local_output_folder, me_at_remote_URL, polyploid_name,
                                   remote_batch_folder):
    if not os.path.exists(local_batch_folder):
        os.makedirs(local_batch_folder)
    run_folders = get_subdirectories(local_output_folder, me_at_remote_URL,
                                     remote_batch_folder, "sp*")
    run_folder_by_polyploid_name = get_run_folders_by_polyploid_name(run_folders)
    print(run_folder_by_polyploid_name)
    remote_path = os.path.join(remote_batch_folder, run_folder_by_polyploid_name[polyploid_name], polyploid_name)
    scp_commands = []
    local_allo_folder = os.path.join(local_batch_folder, polyploid_name)
    if not os.path.exists(local_allo_folder):
        os.makedirs(local_allo_folder)
    remote_to_match = remote_path + "/*final*/*.csv"
    cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_allo_folder]
    print(" ".join(cmd2))
    out_string, error_string = process_wrapper.run_and_wait_with_retry(cmd2, local_allo_folder,
                                                                       "Connection reset by peer", 2, 5)
    time.sleep(5)
    # scp_commands.append(polyploid)
    scp_commands.append(" ".join(cmd2))
    time.sleep(5)
    cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*.xml', '.']
    print(" ".join(cmd2))
    out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)
    print(scp_commands)