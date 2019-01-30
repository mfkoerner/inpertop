# old filtering scheme from get_structures.py


allstructperovtheory = Structure.objects.filter(entry__meta_data__value__contains='perovskite', label='input')
icsd = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input',
    entry__id = F('entry__duplicate_of__id')
    )

# Get all unique icsd entries with 5 atoms and 3 unique species in spacegroup 221
un53icsd221 = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input',
    natoms = 5,
    ntypes = 3,
    entry__id = F('entry__duplicate_of__id'),
    spacegroup = 221
    )
# Knock out all with 3 Oxygen
nonox = [ i for i in un53icsd221 if not 'O3' in i.__str__() ]
# Knock out all with 3 Fluorine
nonoxF = [ i for i in nonox if not 'F3' in i.__str__() ]


nonoxFstr = [i.__str__() for i in nonoxF]

#inverse perovskite, minimum requirement
# invper_min = [i for i in nonoxF if is_inv(i, 'min')]
#inverse perovskite, maximum requirement
invper_max = [i for i in nonoxF if is_inv(i, 'max')]
#inverse perovskite, bond valence requirement
# invper_bv =  [i for i in nonoxF if getv(i, check_inverse = True)]

# filter out F block electrons
invper_nof = [i for i in invper_max if len(PARTIAL_F.intersection(set(i.elements))) == 0]

final_list = invper_nof
# get both types of inverse perovskites (searches for existance of certain wyckoff sites)
# type 1 A: a1; B: b1; X: c3; This is Ram's preferred way
type1 = {i for i in final_list if WSITE_C in [j.wyckoff for j in i.sites]}
# type 2 A: b1; B: a1; X: d3
type2 = {i for i in final_list if WSITE_D in [j.wyckoff for j in i.sites]}
# broken because of bug in qmpy.materials.structures.py definition of translate line 1462
type2_transformed = {struct.recenter(struct[[i.wyckoff.symbol for i in struct.sites].index(u'b')],
 in_place = False, middle = False) for struct in type2}


