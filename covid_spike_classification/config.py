"""Manage configuration options."""

class CSCConfig:
    __slots__ = (
        'datadir',
        'reference',
        'outfile',
        'quiet',
        'show_unexpected',
        'debug',
    )

    def __init__(self, **kwargs):
        for attr in self.__slots__:
            setattr(self, attr, kwargs[attr])

    @classmethod
    def from_args(cls, namespace):
        kwargs = {}
        for attr in cls.__slots__:
            if hasattr(namespace, attr):
                kwargs[attr] = getattr(namespace, attr)

        return cls(**kwargs)