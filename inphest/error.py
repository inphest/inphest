

class InphestException(Exception):
    pass

class FailedSimulationException(InphestException):
    pass

class PreTerminationFailedSimulationException(FailedSimulationException):
    pass

class TotalExtinctionException(PreTerminationFailedSimulationException):
    pass

class PostTerminationFailedSimulationException(FailedSimulationException):
    pass

class SummaryStatisticsCalculationFailure(PostTerminationFailedSimulationException):
    pass

class InsufficientLineagesGenerated(SummaryStatisticsCalculationFailure):
    pass

class IncompleteStateSpaceOccupancyException(SummaryStatisticsCalculationFailure):
    pass

class IncompleteAreaOccupancyException(IncompleteStateSpaceOccupancyException):
    pass

class IncompleteHostOccupancyException(IncompleteStateSpaceOccupancyException):
    pass

