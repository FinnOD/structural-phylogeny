from pathlib import Path

import docker
from docker.types import DeviceRequest


def run_esm_dropouts(
    share_dir: Path, fasta_in: Path, dropout: list[float] | None = None, max_tokens_per_batch: int = 800
) -> Path:

    container = None
    try:
        share_dir_docker = "/home/appuser/bus"
        client = docker.from_env()
        container = client.containers.run(
            "finnod/structphy-esmdropouts-openfold",
            name="structphy-exec",
            detach=True,
            tty=True,
            device_requests=[DeviceRequest(count=-1, capabilities=[["gpu"]])],
            volumes={str(share_dir.resolve()): {"bind": share_dir_docker, "mode": "rw"}},
            user="appuser",
        )

        # Test
        # exec_command = container.exec_run('nvidia-smi')
        # print(exec_command.output.decode())

        # Docker local files
        pybin = "/opt/conda/envs/esmfold/bin/python3"
        infer = "esm-dropouts/scripts/esmfold_inference.py"
        input_fa_full = f"{share_dir_docker}/{str(fasta_in.name)}"
        output_dir_docker = share_dir_docker

        # Make the output dir
        # RAISE if not empty?
        container.exec_run(f"mkdir -p {output_dir_docker}")

        # Run the inference
        inference_command = (
            f"{pybin} {infer} -i {input_fa_full} -o {output_dir_docker} --max-tokens-per-batch {max_tokens_per_batch}"
        )
        if dropout:
            dropout_str = ",".join([f"{x:.4f}" for x in dropout])
            inference_command += f" --dropout {dropout_str}"

        exec_command_gen = container.exec_run(inference_command, stream=True)
        for output in exec_command_gen.output:
            if output:
                print(output.decode().strip())

    finally:
        if container:
            container.stop()
            container.remove()

    return share_dir
