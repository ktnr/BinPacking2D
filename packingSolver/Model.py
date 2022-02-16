
class Item:
    def __init__(self, id, dx, dy):
        self.Id = id
        self.Dx = dx
        self.Dy = dy
        self.Weight = dx * dy
        
    def __eq__(self, other):
        return self.Dx == other.Dx and self.Dy == other.Dy

    def __lt__(self, other):
        return ((self.Weight, self.Dx) < (other.Weight, other.Dx))
        

class Bin:
    def __init__(self, dx, dy):
        self.Dx = dx
        self.Dy = dy
        self.WeightLimit = dx * dy