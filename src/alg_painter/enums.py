import logging
from enum import Enum


class AssemblyMode(Enum):
    AUTO = "auto"
    FINAL = "final"
    DRAFT = "draft"


class GraphType(Enum):
    TILES = "tiles"
    BARS = "bars"
    COMPACT = "compact"
    ALL = "all"


class LogLevel(Enum):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL
