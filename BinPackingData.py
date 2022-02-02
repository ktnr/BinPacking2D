import json

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from Model import Item, Bin

""" from https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html """
def ReadExampleData():
    data = {'bin':{'h':60,'w':40},
            'cat1':{'w': 7,'h':12,'items':10},
            'cat2':{'w': 9,'h': 3,'items':10},
            'cat3':{'w': 5,'h':14,'items':10},
            'cat4':{'w':13,'h': 9,'items':10},
            'cat5':{'w': 6,'h': 8,'items': 5},
            'cat6':{'w':20,'h': 5,'items': 5}}      

    #
    # extract data for easier access
    #

    # bin width and height
    H = data['bin']['h']
    W = data['bin']['w']

    # h,w,cat for each item
    h = [data[cat]['h'] for cat in data if cat!='bin' for i in range(data[cat]['items'])]
    w = [data[cat]['w'] for cat in data if cat!='bin' for i in range(data[cat]['items'])]
    cat = [cat for cat in data if cat!='bin' for i in range(data[cat]['items'])]
    n = len(h)  # number of items
    m = 10  # max number of bins

    return h, w, H, W, m

def ReadBenchmarkData(id):
    data = None
    with open(f'data/CLASS/{id}.json', 'r') as myfile:
        data = myfile.read()

    # parse file
    obj = json.loads(data)

    bin = obj['Objects'][0]
    W = int(bin['Length'])
    H = int(bin['Height'])

    newItems = []
    itemId = 0
    items = obj['Items']
    for item in items:
        for _ in range(0, int(item['Demand'])):
            width = int(item['Length'])
            height = int(item['Height'])

            newItems.append(Item(itemId, width, height))
            itemId += 1

    return newItems, H, W

def ExtractDataForPlot(x, y, items, W, H):
    rectangles = []
    for i, item in enumerate(items):
        solutionX1 = x[i]
        solutionY1 = y[i]
        
        solutionX2 = x[i] + item.Dx
        solutionY2 = y[i] + item.Dy

        color = matplotlib.colors.to_hex([ 1.0 - item.Dx / W, 1.0 - item.Dy / H, 1.0 ])

        rect = [(solutionX1, solutionY1), (solutionX2, solutionY2), color]

        rectangles.append(rect)

    return rectangles

def PlotSolution(maxW, H, rectangles):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #rect1 = matplotlib.patches.Rectangle((0,0), maxW, H, color='c')

    #ax.add_patch(rect1)
    #ax.grid()
    #plt.xlim([0, 20])
    #plt.ylim([0, 6])
    ax.set_xlim((0, maxW))
    ax.set_ylim((0, H))
    ax.set_aspect('equal')

    for i, rect in enumerate(rectangles):
        bottomLeft = rect[0]
        topRight = rect[1]

        x1, y1 = bottomLeft[0], bottomLeft[1]
        x2, y2 = topRight[0], topRight[1]

        color = rect[2]

        rectPlot = matplotlib.patches.Rectangle((x1,y1), x2 - x1, y2 - y1, color=color)
        ax.add_patch(rectPlot)

        ax.annotate(i, (x1 + ((x2 - x1) / 2.0), y1 + ((y2 - y1) / 2.0)), color='b', weight='bold', 
                fontsize=6, ha='center', va='center')

    plt.show()