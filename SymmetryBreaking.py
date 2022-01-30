import math

class SymmetryBreaking:
    @staticmethod
    def DetermineMaximumItemIndexDx(items):
        maximumItemIndex = 0
        maximumItemDx = items[0].Dx
        for i, item in enumerate(items):
            if item.Dx > maximumItemDx:
                maximumItemIndex = i
                maximumItemDx = item.Dx

        return maximumItemIndex

    @staticmethod
    def ReducedDomainX(binDx, item):
        return math.floor((binDx - item.Dx) / 2.0)

    @staticmethod
    def ReducedDomainY(binDy, item):
        return math.floor((binDy - item.Dy) / 2.0)