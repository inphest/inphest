

class InphestException(Exception):
    pass

class FailedSimulationException(InphestException):
    pass

class InsufficientFocalAreaLineagesSimulationException(FailedSimulationException):
    pass

class IncompleteHostOccupancyException(FailedSimulationException):
    pass

class TotalExtinctionException(FailedSimulationException):
    pass

