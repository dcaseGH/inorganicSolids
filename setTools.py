def sameElementByAttributes(obj1, obj2, attributesList):
    ''' Perhaps change this if ever need to test numerics -- only rigorous for strings'''
    return all([obj1.__dict__[a] == obj2.__dict__[a] for a in attributesList])

def subsetByAttributes(objectsList1, objectsList2, attributesList):
    ''' objectsList1 \subseteq objectsList2  '''
    return all([any([sameElementByAttributes(ol1, ol2, attributesList) for ol2 in objectsList2]) for ol1 in objectsList1])

def sameSetByAttributes(objectsList1, objectsList2, attributesList):
    '''   '''
    return len(objectsList1) == len(objectsList2) and subsetByAttributes(objectsList1, objectsList2, attributesList)
