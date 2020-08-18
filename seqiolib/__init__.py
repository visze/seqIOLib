"""Top-level package for seqIOLib."""

__author__ = """Max Schubach"""
__email__ = 'max.schubach@bihealth.de'
__version__ = '0.2.0'

__all__ = ["utils"]

from .sequence import Sequence, Variant, Position, Interval, Orientation, Coordinates, VariantType
from .encoder import Encoder, SequenceCharacterError
