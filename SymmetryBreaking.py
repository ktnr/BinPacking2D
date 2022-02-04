import math

from Model import Item, Bin
from PlacementPoints import PlacementPointGenerator

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

    
    @staticmethod
    def FixIncompatibleItems(incompatibleItems, numberOfItems):
        # This is a similar logic as in section 6.2 in Cote, Haouari, Iori (2019). 
        fixItemToBin = [False] * numberOfItems
        fixItemToBin[0] = True

        if incompatibleItems == None or len(incompatibleItems) == 0:
            return fixItemToBin

        for i in range(1, numberOfItems):
            isIncompatible = True
            for j in range(0, i):
                if frozenset((i, j)) not in incompatibleItems:
                    isIncompatible = False
                    return fixItemToBin
            
            if isIncompatible:
                fixItemToBin[i] = True

        return fixItemToBin

    """ For reference, see section 6.1 from Cote, Iori (2018). """
    @staticmethod
    def CreateReducedBinDomains(incompatibleItems, numberOfItems, numberOfBins, fixItemToBin):
        fullDomains = []
        fullDomains.append([0])
        if incompatibleItems == None or len(incompatibleItems) == 0:
            for i in range(1, numberOfItems):
                boundedNumberOfBins = i if i < numberOfBins else numberOfBins - 1
                fullDomains.append([j for j in range(boundedNumberOfBins + 1)])

            return fullDomains

        
        for i in range(1, numberOfItems):
            boundedNumberOfBins = i if i < numberOfBins else numberOfBins - 1
            fullDomains.append([boundedNumberOfBins])

            if fixItemToBin[i]:
                continue

            for j in range(0, numberOfItems):
                if j >= boundedNumberOfBins:
                    break

                if fixItemToBin[j] and frozenset((i, j)) in incompatibleItems:
                    continue
                else:
                    fullDomains[i].append(j)

        return fullDomains

    def CreateBinDependentNormalPatterns(incompatibleItems, fixItemToBin, items, numberOfBins, binDx, binDy):
        bin = Bin(binDx, binDy)
        numberOfItems = len(items)
        
        itemSpecificNormalPatternsX = []
        itemSpecificNormalPatternsY = []

        # Determine placement points for specific items depending on the combination with every other compatible item.
        for i, itemI in enumerate(items):
            itemSubset = []
            for j, itemJ in enumerate(items):
                if (fixItemToBin[i] and fixItemToBin[j]) or i == j or frozenset((i, j)) in incompatibleItems:
                    continue

                itemSubset.append(itemJ)

            normalPatternsX, normalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy)
            itemSpecificNormalPatternsX.append(normalPatternsX)
            itemSpecificNormalPatternsY.append(normalPatternsY)

        # Determine placement points for items in specific bins and all other compatible items.
        itemBinNormalPatternsX = []
        for i, itemI in enumerate(items):
            itemBinNormalPatternsX.append([])

            if fixItemToBin[i]:
                itemSubset = []
                for j, itemJ in enumerate(items):
                    if fixItemToBin[j] or j <= i or frozenset((i, j)) in incompatibleItems:
                        continue
                    itemSubset.append(itemJ)

                fixedBinNormalPatternsX, fixedBinNormalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy, i*bin.Dx)
                itemBinNormalPatternsX[i].extend(fixedBinNormalPatternsX)

                continue

            for b in range(numberOfBins):
                if b > i:
                    break

                itemSubset = []
                for j, itemJ in enumerate(items):
                    if i == j or b > j or frozenset((i, j)) in incompatibleItems:
                        continue
                    
                    itemSubset.append(itemJ)

                binSpecificNormalPatternsX, binSpecificNormalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy, b*bin.Dx)
                itemBinNormalPatternsX[i].extend(binSpecificNormalPatternsX)

        return itemSpecificNormalPatternsX, itemSpecificNormalPatternsY, itemBinNormalPatternsX