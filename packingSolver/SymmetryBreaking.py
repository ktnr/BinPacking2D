import math

from Model import Item, Bin
from PlacementPoints import PlacementPointGenerator, PlacementPointStrategy
#from packingSolver.PlacementPoints import PlacementPointGenerator, PlacementPointStrategy

class SymmetryBreaking:
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

        placementPointGenerator = PlacementPointGenerator(items, bin, False)

        # Determine placement points for specific items depending on the combination with every other compatible item.
        for i, itemI in enumerate(items):
            itemSubset = []
            for j, itemJ in enumerate(items):
                if (fixItemToBin[i] and fixItemToBin[j]) or i == j or frozenset((i, j)) in incompatibleItems:
                    continue

                itemSubset.append(itemJ)

            #normalPatternsX, normalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy)
            placementPatternsX, placemenetPatternsY = placementPointGenerator.CreatePlacementPatterns(placementPointStrategy, itemI, itemSubset, bin)
            
            if fixItemToBin[i] and PlacementPointGenerator.IsDomainReductionCompatible(placementPointStrategy):
                reducedDomainX = PlacementPointGenerator.ReducedDomainX(bin.Dx, itemI)
                reducedDomainY = PlacementPointGenerator.ReducedDomainY(bin.Dy, itemI)

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

                fixedBinPlacementPatternsX, fixedBinPlacementPatternsY = placementPointGenerator.CreatePlacementPatterns(placementPointStrategy, itemI, itemSubset, bin, i*bin.Dx)
                if PlacementPointGenerator.IsDomainReductionCompatible(placementPointStrategy):
                    reducedDomainFixedBinX = PlacementPointGenerator.ReducedDomainX(bin.Dx, itemI, i*bin.Dx)
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

                binSpecificPlacementPatternsX, binSpecificPlacementPatternsY = placementPointGenerator.CreatePlacementPatterns(placementPointStrategy, itemI, itemSubset, bin, b*bin.Dx)
                if PlacementPointGenerator.IsMaximalPlaceableItem(itemI, itemSubset, domainReducedItems, b) and PlacementPointGenerator.IsDomainReductionCompatible(placementPointStrategy):
                    binSpecificReducedDomainX = PlacementPointGenerator.ReducedDomainX(bin.Dx, itemI, b*bin.Dx)
                    binSpecificPlacementPatternsX = [p for p in binSpecificPlacementPatternsX if p <= binSpecificReducedDomainX]

                    domainReducedItems[itemI.Id].add(b)

                itemBinPlacementPatternsX[i].extend(binSpecificPlacementPatternsX)

        return itemSpecificPlacementPatternsX, itemSpecificPlacementPatternsY, itemBinPlacementPatternsX