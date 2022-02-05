from enum import IntEnum

class PlacementPointStrategy(IntEnum):
    UnitDiscretization = 0
    NormalPatterns = 1
    MeetInTheMiddlePatterns = 2
    MinimalMeetInTheMiddlePatterns = 3

class PlacementPointGenerator:

    @staticmethod
    def CreatePlacementPoints(placementPointStrategy, item, filteredItems, bin):
        binDx = bin.Dx
        binDy = bin.Dy
        placementPointsX = []
        placementPointsY = []
        if placementPointStrategy == PlacementPointStrategy.UnitDiscretization:
            placementPointsX = range(0, binDx + 1 - item.Dx)
            placementPointsY = range(0, binDy + 1 - item.Dy)
        elif placementPointStrategy == PlacementPointStrategy.NormalPatterns:
            placementPointsX, placementPointsY = PlacementPointGenerator.DetermineNormalPatterns(filteredItems, binDx - item.Dx, binDy - item.Dy)
        elif placementPointStrategy == PlacementPointStrategy.MeetInTheMiddlePatterns:
            placementPointsX, placementPointsY = PlacementPointGenerator.DetermineMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy)
            raise ValueError("Meet-in-the-middle patterns might not be accurate, see #2.")
            raise ValueError("Meet-in-the-middle patterns is incompatible with domain reduction.")
        elif placementPointStrategy == PlacementPointStrategy.MinimalMeetInTheMiddlePatterns:
            placementPointsX, placementPointsY = PlacementPointGenerator.DetermineMinimalMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy)
            raise ValueError("Minimal meet-in-the-middle patterns might not be accurate, see #2.")
            raise ValueError("Meet-in-the-middle patterns is incompatible with domain reduction.")
        else:
            raise ValueError('UnkownPlacementPointStrategy')

        return placementPointsX, placementPointsY

    @staticmethod
    def DetermineNormalPatternsX(items, binDx, offsetX = 0):
        if binDx <= 0:
            return [0]

        X = [0] * (binDx + 1)
        X[0] = 1

        for i, item in enumerate(items):
            for p in range(binDx - item.Dx, -1, -1):
                if X[p] == 1:
                    X[p + item.Dx] = 1

        normalPatternsX = []
        for p in range (binDx, -1, -1):
            if X[p] == 1:
                normalPatternsX.append(offsetX + p)

        return normalPatternsX

    @staticmethod
    def DetermineNormalPatternsY(items, binDy):
        if binDy <= 0:
            return [0]

        Y = [0] * (binDy + 1)
        Y[0] = 1

        for i, item in enumerate(items):
            for p in range(binDy - item.Dy, -1, -1):
                if Y[p] == 1:
                    Y[p + item.Dy] = 1

        normalPatternsY = []
        for p in range (binDy, -1, -1):
            if Y[p] == 1:
                normalPatternsY.append(p)

        return normalPatternsY

    @staticmethod
    def DetermineNormalPatterns(items, binDx, binDy, offsetX = 0):
        normalPatternsX = PlacementPointGenerator.DetermineNormalPatternsX(items, binDx)
        normalPatternsY = PlacementPointGenerator.DetermineNormalPatternsY(items, binDy)

        return normalPatternsX, normalPatternsY

    @staticmethod
    def DetermineMeetInTheMiddlePatternsX(items, itemI, binDx, t, offsetX = 0):
        """
        itemI = items[selectedItemIndex]
        filteredItems = []
        for i, item in enumerate(items):
            if i == selectedItemIndex:
                continue
            filteredItems.append(item)
        """

        meetInTheMiddlePoints = PlacementPointGenerator.DetermineNormalPatternsX(items, min(t - 1, binDx - itemI.Dx), offsetX) # placemenetPointsLeft
        placemenetPointsRightPrime = PlacementPointGenerator.DetermineNormalPatternsX(items, binDx - itemI.Dx - t, offsetX)

        for p in placemenetPointsRightPrime:
            meetInTheMiddlePoints.append(binDx - itemI.Dx - p)

        return meetInTheMiddlePoints

    @staticmethod
    def DetermineMeetInTheMiddlePatternsY(items, itemI, binDy, t):
        """
        itemI = items[selectedItemIndex]
        filteredItems = []
        for i, item in enumerate(items):
            if i == selectedItemIndex:
                continue
            filteredItems.append(item)
        """

        meetInTheMiddlePoints = PlacementPointGenerator.DetermineNormalPatternsY(items, min(t - 1, binDy - itemI.Dy)) # placemenetPointsLeft
        placemenetPointsRightPrime = PlacementPointGenerator.DetermineNormalPatternsY(items, binDy - itemI.Dy - t)

        for p in placemenetPointsRightPrime:
            meetInTheMiddlePoints.append(binDy - itemI.Dy - p)

        return meetInTheMiddlePoints

    @staticmethod
    def DetermineMeetInTheMiddlePatterns(items, itemI, binDx, binDy, offsetX = 0):
        meetInTheMiddlePointsX = set()
        meetInTheMiddlePointsY = set()
        for t in range(1, binDx + 1):
            meetInTheMiddlePoints = PlacementPointGenerator.DetermineMeetInTheMiddlePatternsX(items, itemI, binDx, t, offsetX)
            #meetInTheMiddlePointsX.extend(meetInTheMiddlePoints)
            meetInTheMiddlePointsX.update(meetInTheMiddlePoints)

        for t in range(1, binDy + 1):
            meetInTheMiddlePoints = PlacementPointGenerator.DetermineMeetInTheMiddlePatternsY(items, itemI, binDy, t)
            #meetInTheMiddlePointsY.extend(meetInTheMiddlePoints)
            meetInTheMiddlePointsY.update(meetInTheMiddlePoints)

        return list(meetInTheMiddlePointsX), list(meetInTheMiddlePointsY)
        
    @staticmethod
    def DetermineMinimalMeetInTheMiddlePatternsX(items, itemI, binDx, offsetX = 0):
        meetInTheMiddlePointsLeftX = [0] * (binDx + 1)
        meetInTheMiddlePointsRightX = [0] * (binDx + 1)

        regularNormalPatterns = PlacementPointGenerator.DetermineNormalPatternsX(items, binDx - itemI.Dx, offsetX)
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
        
    @staticmethod
    def DetermineMinimalMeetInTheMiddlePatternsY(items, itemI, binDy):
        meetInTheMiddlePointsLeftY = [0] * (binDy + 1)
        meetInTheMiddlePointsRightY = [0] * (binDy + 1)

        normalPatternsY = PlacementPointGenerator.DetermineNormalPatternsY(items, binDy - itemI.Dy)

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
        
    @staticmethod
    def DetermineMinimalMeetInTheMiddlePatterns(items, itemI, binDx, binDy, offsetX = 0):
        meetInTheMiddlePointsX = PlacementPointGenerator.DetermineMinimalMeetInTheMiddlePatternsX(items, itemI, binDx, offsetX)
        meetInTheMiddlePointsY = PlacementPointGenerator.DetermineMinimalMeetInTheMiddlePatternsY(items, itemI, binDy)

        return meetInTheMiddlePointsX, meetInTheMiddlePointsY