import base64
import subprocess
import sys

import click


@click.command()
@click.argument("name")
def hello(name: str):
    if sys.version_info.minor != 11:
        raise RuntimeError("Wrong python version")

    path = subprocess.run(["which", "conda"], capture_output=True).stdout.decode()
    if "forge" not in path or "anaconda" in path:
        raise RuntimeError("Wrong conda version")

    subprocess.run(["conda", "list"], check=True, capture_output=True)

    print("Your answer is:")
    eval(
        base64.b64decode(
            b"cHJpbnQobGVuKG5hbWUubG93ZXIoKSkgKyBzeXMudmVyc2lvbl9pbmZvLm1pbm9yKQ=="
        )
    )


if __name__ == "__main__":
    hello()
