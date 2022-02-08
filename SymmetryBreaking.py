import math

from Model import Item, Bin
from PlacementPoints import PlacementPointGenerator, PlacementPointStrategy

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

    @staticmethod
    def DetermineFixedItems(incompatibleItems, numberOfItems):
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

    @staticmethod
    def CreateBinDependentPlacementPatterns(incompatibleItems, fixItemToBin, items, numberOfBins, bin, placementPointStrategy = PlacementPointStrategy.NormalPatterns):        
        itemSpecificPlacementPatternsX = []
        itemSpecificPlacementPatternsY = []

        # Determine placement points for specific items depending on the combination with every other compatible item.
        for i, itemI in enumerate(items):
            itemSubset = []
            for j, itemJ in enumerate(items):
                if (fixItemToBin[i] and fixItemToBin[j]) or i == j or frozenset((i, j)) in incompatibleItems:
                    continue

                itemSubset.append(itemJ)

            #normalPatternsX, normalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy)
            placementPatternsX, placemenetPatternsY = PlacementPointGenerator.CreatePlacementPatterns(placementPointStrategy, itemI, itemSubset, bin)
            
            if fixItemToBin[i] and SymmetryBreaking.IsDomainReductionCompatible(placementPointStrategy):
                reducedDomainX = SymmetryBreaking.ReducedDomainX(bin.Dx, itemI)
                reducedDomainY = SymmetryBreaking.ReducedDomainY(bin.Dy, itemI)

                itemSpecificPlacementPatternsX.append([p for p in placementPatternsX if p <= reducedDomainX])
                itemSpecificPlacementPatternsY.append([p for p in placemenetPatternsY if p <= reducedDomainY])
            else:
                itemSpecificPlacementPatternsX.append(placementPatternsX)
                itemSpecificPlacementPatternsY.append(placemenetPatternsY)

        # Determine placement points for items in specific bins and all other compatible items.
        domainReducedItems = [set() for _ in range(len(items))]
        itemBinPlacementPatternsX = []
        for i, itemI in enumerate(items):
            itemBinPlacementPatternsX.append([])

            if fixItemToBin[i]:
                itemSubset = []
                for j, itemJ in enumerate(items):
                    if fixItemToBin[j] or j <= i or frozenset((i, j)) in incompatibleItems:
                        continue
                    itemSubset.append(itemJ)

                fixedBinPlacementPatternsX, fixedBinPlacementPatternsY = PlacementPointGenerator.CreatePlacementPatterns(placementPointStrategy, itemI, itemSubset, bin, i*bin.Dx)
                if SymmetryBreaking.IsDomainReductionCompatible(placementPointStrategy):
                    reducedDomainFixedBinX = SymmetryBreaking.ReducedDomainX(bin.Dx, itemI, i*bin.Dx)
                    fixedBinPlacementPatternsX = [p for p in fixedBinPlacementPatternsX if p <= reducedDomainFixedBinX]

                    domainReducedItems[itemI.Id].add(i)

                itemBinPlacementPatternsX[i].extend(fixedBinPlacementPatternsX)

                continue

            for b in range(numberOfBins):
                if b > i:
                    break

                itemSubset = []
                for j, itemJ in enumerate(items):
                    if i == j or b > j or frozenset((i, j)) in incompatibleItems:
                        continue
                    
                    itemSubset.append(itemJ)

                binSpecificPlacementPatternsX, binSpecificPlacementPatternsY = PlacementPointGenerator.CreatePlacementPatterns(placementPointStrategy, itemI, itemSubset, bin, b*bin.Dx)
                if SymmetryBreaking.IsMaximalPlaceableItem(itemI, itemSubset, domainReducedItems, b) and SymmetryBreaking.IsDomainReductionCompatible(placementPointStrategy):
                    binSpecificReducedDomainX = SymmetryBreaking.ReducedDomainX(bin.Dx, itemI, b*bin.Dx)
                    binSpecificPlacementPatternsX = [p for p in binSpecificPlacementPatternsX if p <= binSpecificReducedDomainX]

                    domainReducedItems[itemI.Id].add(b)

                itemBinPlacementPatternsX[i].extend(binSpecificPlacementPatternsX)

        return itemSpecificPlacementPatternsX, itemSpecificPlacementPatternsY, itemBinPlacementPatternsX