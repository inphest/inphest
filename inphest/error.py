

class InphestException(Exception):
    pass

class FailedSimulationException(InphestException):
    pass

class PostTerminationFailedSimulationException(FailedSimulationException):
    pass

class PreTerminationFailedSimulationException(FailedSimulationException):
    pass

class IncompleteAreaOccupancyException(PostTerminationFailedSimulationException):
    pass

class IncompleteHostOccupancyException(PostTerminationFailedSimulationException):
    pass

class TotalExtinctionException(PreTerminationFailedSimulationException):
    pass

