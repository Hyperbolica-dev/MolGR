import importlib.metadata


try:
    __version__ = importlib.metadata.version("myrepositorytemplate")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
