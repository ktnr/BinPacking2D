from enum import IntEnum
from Model import Item

#import SymmetryBreaking

import math

class PlacementPointStrategy(IntEnum):
    StandardUnitDiscretization = 0
    UnitDiscretization = 1
    NormalPatterns = 2
    MeetInTheMiddlePatterns = 3
    MinimalMeetInTheMiddlePatterns = 4

class PlacementPointGenerator:
    def __init__(self, items, bin, enableDomainReduction = True):
        if enableDomainReduction:
            # Symmetry breaking according to 
            # Soh, T., Inoue, K., Tamura, N., Banbara, M., & Nabeshima, H. (2010). A SAT-based method for solving the two-dimensional strip packing problem. Fundamenta Informaticae, 102(3-4), 467-487.
            self.reducedItemIndex = PlacementPointGenerator.DetermineMaximumItemIndexDx(items)

            self.reducedItem = items[self.reducedItemIndex]
            self.reducedDomainThresholdX = PlacementPointGenerator.ReducedDomainX(bin.Dx, self.reducedItem)
            self.reducedDomainThresholdY = PlacementPointGenerator.ReducedDomainY(bin.Dy, self.reducedItem)
        else:
            self.reducedItemIndex = -1

            self.reducedItem = Item(-1, 0, 0)
            self.reducedDomainThresholdX = bin.Dx
            self.reducedDomainThresholdY = bin.Dy

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
    def ReducedDomainX(binDx, item, offsetX = 0):
        return offsetX + math.floor((binDx - item.Dx) / 2.0)

    @staticmethod
    def ReducedDomainY(binDy, item):
        return math.floor((binDy - item.Dy) / 2.0)

    @staticmethod
    def IsDomainReductionCompatible(placementPointStrategy):
        if placementPointStrategy == PlacementPointStrategy.UnitDiscretization or placementPointStrategy == PlacementPointStrategy.NormalPatterns:
            return True
        
        return False
    
    @staticmethod
    def IsMaximalPlaceableItem(item, itemSubset, domainReducedItems, binId):
        for otherItem in itemSubset:
            if otherItem.Dx > item.Dx or binId in domainReducedItems[otherItem.Id]:
            #if otherItem.Dx > item.Dx or domainReducedItems[otherItem.Id] == True:
                return False

        return True

    def CreatePlacementPatterns(self, placementPointStrategy, item, filteredItems, bin, offsetX = 0):
        # TODO.Logic: bring for-loop over items into this methods and return a list of placement patterns, one entry for each item.
        # Before the loop, check if placementPointStrategy == PlacementPointStrategy.MinimalMeetInTheMiddlePatterns. If so, jump into separate routine and return immediately afterwards. 
        binDx = bin.Dx
        binDy = bin.Dy
        placementPointsX = []
        placementPointsY = []
        if placementPointStrategy == PlacementPointStrategy.StandardUnitDiscretization:
            placementPointsX, placementPointsY = self.GenerateStandardUnitDiscretization(filteredItems, item, bin, offsetX)
        elif placementPointStrategy == PlacementPointStrategy.UnitDiscretization:
            placementPointsX, placementPointsY = self.GenerateUnitDiscretization(filteredItems, item, bin, offsetX)
        elif placementPointStrategy == PlacementPointStrategy.NormalPatterns:
            placementPointsX, placementPointsY = self.GenerateNormalPatterns(filteredItems, item, binDx - item.Dx, binDy - item.Dy, offsetX)
        elif placementPointStrategy == PlacementPointStrategy.MeetInTheMiddlePatterns:
            placementPointsX, placementPointsY = self.GenerateMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy, offsetX)
        elif placementPointStrategy == PlacementPointStrategy.MinimalMeetInTheMiddlePatterns:
            placementPointsX, placementPointsY = self.GenerateMinimalMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy, offsetX)
            raise ValueError("Minimal meet-in-the-middle patterns are not accurate, see #2.")
        else:
            raise ValueError('UnkownPlacementPointStrategy')

        return placementPointsX, placementPointsY

    """ Without symmetry breaking. """
    def GenerateStandardUnitDiscretization(self, filteredItems, item, bin, offsetX = 0):
        placementPointsX = list(range(offsetX, offsetX + bin.Dx - item.Dx + 1))
        placementPointsY = list(range(0, bin.Dy - item.Dy + 1))

        return placementPointsX, placementPointsY

    """ With symmetry breaking. """
    def GenerateUnitDiscretization(self, filteredItems, item, bin, offsetX = 0):
        placementEndX = offsetX + bin.Dx - item.Dx
        placementEndY =  bin.Dy - item.Dy
        if item.Id == self.reducedItem.Id:
            placementEndX = min(placementEndX, offsetX + self.reducedDomainThresholdX)
            placementEndY = min(placementEndY, self.reducedDomainThresholdY)

        placementPointsX = list(range(offsetX, placementEndX + 1))
        placementPointsY = list(range(0, placementEndY + 1))

        return placementPointsX, placementPointsY

    """ 
    Domain reduction is feasible for normal patterns. Any feasible solution where the reduced item is not placed within its reduced
    domain can be transformed into a solution where it is placed in its reduced domain by mirroring the solution.
    """
    def DetermineNormalPatternsX(self, items, item, binDx, offsetX = 0):
        if binDx <= 0:
            return [offsetX]

        X = [0] * (binDx + 1)
        X[0] = 1

        for i, item in enumerate(items):
            for p in range(binDx - item.Dx, -1, -1):
                if X[p] == 1:
                    X[p + item.Dx] = 1

        normalPatternsX = []
        decrementStart = binDx
        if item.Id == self.reducedItem.Id:
            decrementStart = self.reducedDomainThresholdX

        for p in range (decrementStart, -1, -1):
            if X[p] == 1:
                normalPatternsX.append(offsetX + p)

        return normalPatternsX

    def DetermineNormalPatternsY(self, items, item, binDy):
        if binDy <= 0:
            return [0]

        Y = [0] * (binDy + 1)
        Y[0] = 1

        for i, item in enumerate(items):
            for p in range(binDy - item.Dy, -1, -1):
                if Y[p] == 1:
                    Y[p + item.Dy] = 1

        normalPatternsY = []
        decrementStart = binDy
        if item.Id == self.reducedItem.Id:
            decrementStart = self.reducedDomainThresholdY

        for p in range (decrementStart, -1, -1):
            if Y[p] == 1:
                normalPatternsY.append(p)

        return normalPatternsY

    """
    # In this variant, the fact that the item with a reduced domain can only be placed in a certain range is also considered
    # when generating placement points for other items, not only for the reduced item itself.
    # This might preserve optimality/feasibility but requires proof. What is special with this variant is that the order in 
    # which the placement points are generated matters, i.e., the order of the items list, which was not the case before. 
    # Consider a bin of dimension 10 and the item dimensions 1, 4, and 5. Reduced normal patterns are:
    # 5: 0, 1, 4, 5 --> reduced domain = floor((10 - 5) / 2) = 2 --> 0, 1
    # 4: 0, 1, 5, 6
    # 1 (5, 4): 0, 4, 5, 9
    # 1 (4, 5): 0, 4, 5 (9 is not generated because the reduced item cannot be placed on 4)
    # With order (5, 4), 
    # - the order 5, 4, 1 is possible so that item 1 can be placed at 9,
    # - the order 4, 5, 1 is not possible because item 5 is outside its reduced domain.
    # - the order 1, 5, 4 is possible and the mirror solution to 4, 5, 1
    # With order (4, 5), 
    # - the order 5, 4, 1, so that item 1 cannot be placed at 9, is not possible. 
    # - the order 4, 5, 1 is not possible
    # - the order 1, 5, 4 is possible and its mirror solution 4, 5, 1 hast item 1 is placed at 9.
    # Does this exclude feasible solutions?
    def DetermineNormalPatternsX(self, items, item, binDx, offsetX = 0):
        if binDx <= 0:
            return [offsetX]

        X = [0] * (binDx + 1)
        X[0] = 1

        for i, item in enumerate(items):
            itemDecrementStart = binDx - item.Dx
            if item.Id == self.reducedItem.Id:
                itemDecrementStart = min(itemDecrementStart, self.reducedDomainThresholdX)

            for p in range(itemDecrementStart, -1, -1):
                if X[p] == 1:
                    X[p + item.Dx] = 1

        normalPatternsX = []
        decrementStart = binDx
        if item.Id == self.reducedItem.Id:
            decrementStart = self.reducedDomainThresholdX

        for p in range (decrementStart, -1, -1):
            if X[p] == 1:
                normalPatternsX.append(offsetX + p)

        return normalPatternsX

    def DetermineNormalPatternsY(self, items, item, binDy):
        if binDy <= 0:
            return [0]

        Y = [0] * (binDy + 1)
        Y[0] = 1

        for i, item in enumerate(items):
            itemDecrementStart = binDy - item.Dy
            if item.Id == self.reducedItem.Id:
                itemDecrementStart = min(itemDecrementStart, self.reducedDomainThresholdY)

            for p in range(itemDecrementStart, -1, -1):
                if Y[p] == 1:
                    Y[p + item.Dy] = 1

        normalPatternsY = []
        decrementStart = binDy
        if item.Id == self.reducedItem.Id:
            decrementStart = self.reducedDomainThresholdY

        for p in range (decrementStart, -1, -1):
            if Y[p] == 1:
                normalPatternsY.append(p)

        return normalPatternsY
    """

    def GenerateNormalPatterns(self, items, item, binDx, binDy, offsetX = 0):
        normalPatternsX = self.DetermineNormalPatternsX(items, item, binDx, offsetX)
        normalPatternsY = self.DetermineNormalPatternsY(items, item, binDy)

        return normalPatternsX, normalPatternsY

    """ 
    Domain reduction on normal patterns is compatible with meet-in-the-middle patterns. The proof should be similar to that of Proposition 5
    in Cote and Iori (2018): Meet-in-the-middle principle, which more or less states that one item can be placed in one corner of the container.
    """
    def DetermineMeetInTheMiddlePatternsX(self, items, itemI, binDx, t, offsetX = 0):
        """
        itemI = items[selectedItemIndex]
        filteredItems = []
        for i, item in enumerate(items):
            if i == selectedItemIndex:
                continue
            filteredItems.append(item)
        """

        meetInTheMiddlePoints = self.DetermineNormalPatternsX(items, itemI, min(t - 1, binDx - itemI.Dx), offsetX) # placemenetPointsLeft
        placemenetPointsRightPrime = self.DetermineNormalPatternsX(items, itemI, binDx - itemI.Dx - t, offsetX)

        for p in placemenetPointsRightPrime:
            meetInTheMiddlePoints.append(offsetX + binDx - itemI.Dx - p)

        return meetInTheMiddlePoints

    def DetermineMeetInTheMiddlePatternsY(self, items, itemI, binDy, t):
        """
        itemI = items[selectedItemIndex]
        filteredItems = []
        for i, item in enumerate(items):
            if i == selectedItemIndex:
                continue
            filteredItems.append(item)
        """

        meetInTheMiddlePoints = self.DetermineNormalPatternsY(items, itemI, min(t - 1, binDy - itemI.Dy)) # placemenetPointsLeft
        placemenetPointsRightPrime = self.DetermineNormalPatternsY(items, itemI, binDy - itemI.Dy - t)

        for p in placemenetPointsRightPrime:
            meetInTheMiddlePoints.append(binDy - itemI.Dy - p)

        return meetInTheMiddlePoints

    def GenerateMeetInTheMiddlePatterns(self, items, itemI, binDx, binDy, offsetX = 0):
        meetInTheMiddlePointsX = set()
        meetInTheMiddlePointsY = set()
        for t in range(1, binDx + 1):
            meetInTheMiddlePoints = self.DetermineMeetInTheMiddlePatternsX(items, itemI, binDx, t, offsetX)
            #meetInTheMiddlePointsX.extend(meetInTheMiddlePoints)
            meetInTheMiddlePointsX.update(meetInTheMiddlePoints)

        for t in range(1, binDy + 1):
            meetInTheMiddlePoints = self.DetermineMeetInTheMiddlePatternsY(items, itemI, binDy, t)
            #meetInTheMiddlePointsY.extend(meetInTheMiddlePoints)
            meetInTheMiddlePointsY.update(meetInTheMiddlePoints)

        return list(meetInTheMiddlePointsX), list(meetInTheMiddlePointsY)
        
    def DetermineMinimalMeetInTheMiddlePatternsX(self, items, itemI, binDx, offsetX = 0):
        meetInTheMiddlePointsLeftX = [0] * (binDx + 1)
        meetInTheMiddlePointsRightX = [0] * (binDx + 1)

        regularNormalPatterns = self.DetermineNormalPatternsX(items, itemI, binDx - itemI.Dx, offsetX)
        #normalPatternsX = PlacementPointGenerator.DetermineNormalPatternsX(items, binDx - itemI.Dx, offsetX)

        #regularNormalPatterns = set(normalPatternsX)
        for p in regularNormalPatterns:
            meetInTheMiddlePointsLeftX[p] = 1 # meetInTheMiddlePointsLeftX[p] + 1
            meetInTheMiddlePointsRightX[binDx - itemI.Dx - p] = 1

        for p in range(1, binDx + 1):
            meetInTheMiddlePointsLeftX[p] = meetInTheMiddlePointsLeftX[p] + meetInTheMiddlePointsLeftX[p - 1]
            meetInTheMiddlePointsRightX[binDx - p] = meetInTheMiddlePointsRightX[binDx - p] + meetInTheMiddlePointsRightX[binDx - (p - 1)]

        tMin = 1
        minX = meetInTheMiddlePointsLeftX[0] + meetInTheMiddlePointsRightX[1]

        for p in range(2, binDx + 1):
            if meetInTheMiddlePointsLeftX[p - 1] + meetInTheMiddlePointsRightX[p] < minX:
                minX = meetInTheMiddlePointsLeftX[p - 1] + meetInTheMiddlePointsRightX[p]
                tMin = p

        meetInTheMiddlePointsX = []
        for p in regularNormalPatterns:
            if p < tMin:
                meetInTheMiddlePointsX.append(p)

            if binDx - itemI.Dx - p >= tMin:
                meetInTheMiddlePointsX.append(binDx - itemI.Dx - p)

        return meetInTheMiddlePointsX
        
    def DetermineMinimalMeetInTheMiddlePatternsY(self, items, itemI, binDy):
        meetInTheMiddlePointsLeftY = [0] * (binDy + 1)
        meetInTheMiddlePointsRightY = [0] * (binDy + 1)

        normalPatternsY = self.DetermineNormalPatternsY(items, itemI, binDy - itemI.Dy)

        regularNormalPatterns = set(normalPatternsY)
        for p in regularNormalPatterns:
            meetInTheMiddlePointsLeftY[p] = 1 # meetInTheMiddlePointsLeftX[p] + 1
            meetInTheMiddlePointsRightY[binDy - itemI.Dy - p] = 1

        for p in range(1, binDy + 1):
            meetInTheMiddlePointsLeftY[p] = meetInTheMiddlePointsLeftY[p] + meetInTheMiddlePointsLeftY[p - 1]
            meetInTheMiddlePointsRightY[binDy - p] = meetInTheMiddlePointsRightY[binDy - p] + meetInTheMiddlePointsRightY[binDy - (p - 1)]

        tMin = 1
        minY = meetInTheMiddlePointsLeftY[0] + meetInTheMiddlePointsRightY[1]

        for p in range(2, binDy + 1):
            if meetInTheMiddlePointsLeftY[p - 1] + meetInTheMiddlePointsRightY[p] < minY:
                minY = meetInTheMiddlePointsLeftY[p - 1] + meetInTheMiddlePointsRightY[p]
                tMin = p

        meetInTheMiddlePointsY = []
        for p in regularNormalPatterns:
            if p < tMin:
                meetInTheMiddlePointsY.append(p)

            if binDy - itemI.Dy - p >= tMin:
                meetInTheMiddlePointsY.append(binDy - itemI.Dy - p)

        return meetInTheMiddlePointsY
        
    def GenerateMinimalMeetInTheMiddlePatterns(self, items, itemI, binDx, binDy, offsetX = 0):
        meetInTheMiddlePointsX = self.DetermineMinimalMeetInTheMiddlePatternsX(items, itemI, binDx, offsetX)
        meetInTheMiddlePointsY = self.DetermineMinimalMeetInTheMiddlePatternsY(items, itemI, binDy)

        return meetInTheMiddlePointsX, meetInTheMiddlePointsY