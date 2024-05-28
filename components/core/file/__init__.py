""" Provides helpers for working with files

Newly created files should be upserted with the upsert_file function,
this will ensure the file is accessible to any worker that needs to
access it.
"""

import os
import typing
import tempfile
import requests
import contextlib

from ufs.access.url import register_proto_handler


def Storage():
    import json
    from ufs.spec import UFS

    if "UFS_STORAGE" in os.environ:
        return UFS.from_dict(**json.loads(os.environ["UFS_STORAGE"]))
    else:
        import pathlib
        from ufs.impl.local import Local
        from ufs.impl.prefix import Prefix

        storagePath = (
            pathlib.Path(os.environ.get("PYTHON_ROOT", pathlib.Path.cwd()))
            / "data"
            / "ufs"
        )
        storagePath.mkdir(parents=True, exist_ok=True)
        return Prefix(Local(), storagePath)


@register_proto_handler("storage")
def proto_storage(url):
    from ufs.impl.prefix import Prefix

    return Prefix(Storage(), url["path"])


class File(typing.TypedDict):
    url: str
    filename: str
    description: typing.Optional[str]
    size: int
    sha256: typing.Optional[str]


class TemporaryFile(dict):
    def __init__(self, suffix="", description=None):
        super().__init__(
            url="",
            filename="",
            description=description,
            size=0,
        )
        self.file = tempfile.mktemp(suffix=suffix)
        self["filename"] = self.file

    def __repr__(self) -> str:
        return f"TemporaryFile({repr(self)})"

    def unlink(self):
        os.unlink(self.file)


@contextlib.contextmanager
def upsert_file(suffix="", description=None):
    """Usage:
    from components.core.file import upsert_file
    with upsert_file('.csv') as f:
      # or, for example: df.to_csv(f.file)
    # f.url contains the url to the file uploaded
    print(f.url)
    """
    # NOTE: this context manager will remove the temporary file when done
    tmp = TemporaryFile(suffix, description=description)
    try:
        # we give you a temporary file to write to
        yield tmp
        # we upload the file to the main server to get a persistent url
        req = requests.post(
            f"{os.environ['PUBLIC_URL']}/api/v1/components/core/file/upload",
            headers={"Authorization": f"Token {os.environ['NEXTAUTH_SECRET']}"},
            files=[("file", open(tmp.file, "rb"))],
        )
        assert (
            req.status_code >= 200 and req.status_code < 300
        ), f"Error ({req.status_code}): {req.json()}"
        res = req.json()
        # we return register url to the file
        tmp.update(res["file"][0])
    finally:
        # we remove the temporary file
        tmp.unlink()


@contextlib.contextmanager
def file_as_stream(file: File, *args, **kwargs) -> str:
    from ufs.access.url import open_from_url

    with open_from_url(file, *args, **kwargs) as fr:
        yield fr


@contextlib.contextmanager
def file_as_path(file: File, *args, **kwargs) -> str:
    from ufs.access.url import ufs_file_from_url
    from ufs.access.ffuse import ffuse_mount

    ufs, filename = ufs_file_from_url(file["url"], filename=file["filename"])
    with ffuse_mount(ufs, readonly=True) as mount_dir:
        yield str(mount_dir / filename)


def file_read_stream(file: File, chunk_size=8192):
    with file_as_stream(file, "rb") as fr:
        while True:
            buf = fr.read(chunk_size)
            if not buf:
                break
            yield buf


def file_move(origFile: File, newFile: File):
    from ufs.access.url import ufs_file_from_url
    from ufs.access.shutil import movefile

    orig_ufs, orig_file = ufs_file_from_url(origFile["url"])
    new_ufs, new_file = ufs_file_from_url(newFile["url"])
    with orig_ufs:
        with new_ufs:
            movefile(orig_ufs, orig_file, new_ufs, new_file)


def one_or_none(it):
    try:
        one = next(it)
    except StopIteration:
        return None
    try:
        next(it)
    except StopIteration:
        return one
    else:
        return None


def file_stats_from_url(url: str, chunk_size=8192):
    import hashlib

    h = hashlib.sha256()
    s = 0
    with file_as_stream(dict(url=url), "rb") as fr:
        while True:
            buf = fr.read(chunk_size)
            if not buf:
                break
            h.update(buf)
            s += len(buf)
    return dict(
        size=s,
        sha256=h.hexdigest(),
    )


def file_from_url(url: str):
    from ufs.utils.url import parse_url, parse_netloc

    url_parsed = parse_url(url)
    if url_parsed["proto"] == "drs":
        import requests

        netloc_parsed = parse_netloc(url_parsed)
        scheme = "http" if netloc_parsed["port"] else "https"
        req = requests.get(
            f"{scheme}://{netloc_parsed['netloc']}/ga4gh/drs/v1/objects{netloc_parsed['path']}"
        )
        assert req.ok
        file_info = req.json()
        filename = file_info.get("name") or url_parsed["path"].rpartition("/")[-1]
        sha256 = one_or_none(
            checksum["checksum"]
            for checksum in file_info.get("checksums", [])
            if checksum.get("type") == "sha-256"
        )
        if sha256:
            size = file_info["size"]
            return dict(url=url, filename=filename, size=size, sha256=sha256)
        else:
            return dict(url=url, filename=filename, **file_stats_from_url(url))
    else:
        return dict(
            url=url,
            filename=url_parsed["path"].rpartition("/")[-1],
            **file_stats_from_url(url),
        )
