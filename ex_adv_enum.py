from enum import Enum

class NoValue(Enum):
    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)

class AutoNumber(NoValue):
    def __new__(cls, init=None):
        _value: int
        if isinstance(init, int):
            _value = init
        elif isinstance(init, str):
            _value = cls[init].value
        else:
            _value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = _value
        return obj

class Color(AutoNumber):
    RED = ()
    GREEN = 5
    BLUE = ()
    MARS = 'RED'
    APPLE = 'GREEN'


print(Color.RED.value)
print(Color.GREEN.value)
print(Color.BLUE.value)
print(Color.MARS.value)
print(Color.APPLE.value)
