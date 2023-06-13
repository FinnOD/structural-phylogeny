from pathlib import Path
from typing import List
import docker

def run_esm_dropouts(share_dir: Path, fasta_in: Path, dropout: List[float] = None, max_tokens_per_batch: int =800) -> List[Path]: #returns list[pdb_file, pdb_file] 
    try:
        share_dir_docker = '/mnt/bus'

        client = docker.from_env()
        print(client)
        container = client.containers.run(
            'finnod/structphy-esmdropouts-openfold',
            name='structphy-exec',
            detach=True,
            tty=True,
            device_requests=[docker.types.DeviceRequest(count=-1, capabilities=[['gpu']])],
            volumes={share_dir: {'bind': share_dir_docker, 'mode': 'rw'}},
        )
        print('container', container)
        # Test
        exec_command = container.exec_run('nvidia-smi')
        print(exec_command.output.decode())

        
        # Docker local files
        pybin = '/opt/conda/envs/esmfold/bin/python3'
        infer = 'esm-dropouts/scripts/esmfold_inference.py'
        input_fa_full = f'{share_dir_docker}/{str(fasta_in)}'
        output_dir_docker = f'{share_dir_docker}/pdbout/'
        
        # Make the output dir
        exec_command = container.exec_run(f'mkdir -p {output_dir_docker}')
        
        # Run the inference
        inference_command = f'{pybin} {infer} -i {input_fa_full} -o {output_dir_docker} --max-tokens-per-batch {max_tokens_per_batch} '
        if dropout:
            dropout_str = ",".join([f'{x:.4f}' for x in dropout])
            inference_command += f'--dropout {dropout_str}'
        
        # print(inference command:
        # print(inference_command)
        # exec_command_gen = container.exec_run(inference_command, stream=True)
        # for output in exec_command_gen.output:
        #     if output:
        #         print(output.decode().strip())

        exec_command = container.exec_run(f'ls {output_dir_docker}')
        out_ls = exec_command.output.decode().strip().split('\n')
        return [share_dir+'pdbout/'+file for file in out_ls]

        # Close and delete the container
    except Exception as e:
        raise e
    finally:
        container.stop()
        container.remove()
        

a = run_esm_dropouts('/home/ubuntu/bus/', 'single.fa', dropout=[0.99]*36)
print(a)