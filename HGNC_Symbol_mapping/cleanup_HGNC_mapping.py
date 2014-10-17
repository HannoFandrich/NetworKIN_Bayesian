from Lindinglab.idmapper import *



manual_map = {
    "ENSP00000350696": "ZNF587B",   # Uniprot
    "ENSP00000299353": "C10orf32",   # Uniprot
    "ENSP00000341497": "ZNF177",   # Uniprot
    "ENSP00000382819": "DXO",   # Uniprot
    "ENSP00000395262": "DXO",   # Uniprot
    "ENSP00000410404": "DXO",   # Uniprot
    "ENSP00000226218": "VTN",   # Uniprot
    "ENSP00000391123": "DXO",   # Uniprot
    "ENSP00000391349": "DXO",   # Uniprot
    "ENSP00000383645": "KIR2DL5A",   # Uniprot
    "ENSP00000301407": "CGB1",   # Uniprot
    "ENSP00000316021": "PLSCR3",   # Uniprot
    }


idmap_Biomart_Ens75 = ReadIdMapSimple("Biomart_map_Ensembl75.txt", 1, 3, offset = 1)
idmap_Biomart_Ens68 = ReadIdMapSimple("Biomart_map_Ensembl68.txt", 1, 3, offset = 1)


idmap_Biomart_Ensp_HGNC = CIdMap()

MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl75.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl68.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl67.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl66.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl65.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl64.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl63.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl62.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl59.txt", 1, 3, offset = 1))
MergeMap(idmap_Biomart_Ensp_HGNC, ReadIdMapSimple("Biomart_map_Ensembl54.txt", 1, 3, offset = 1))

idmap_Biomart_Ensp_Ensg = CIdMap()

MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl75.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl68.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl67.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl66.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl65.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl64.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl63.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl62.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl59.txt", 1, 2, offset = 1))
MergeMap(idmap_Biomart_Ensp_Ensg, ReadIdMapSimple("Biomart_map_Ensembl54.txt", 1, 2, offset = 1))


idmap_HGNC_Ensg_HGNC = ReadIdMapSimple("idmap_HGNC.txt", 8, 2, offset = 1)
MergeMap(idmap_HGNC_Ensg_HGNC, ReadIdMapSimple("idmap_HGNC.txt", 10, 2, offset = 1))


idmap_Uniprot_StringID_UniprotAC = ReadIdMapSimple("idmap_STRING_UniprotAC_by_Uniprot.txt")
idmap_

approved_symbol_names = set(ReadIdMapSimple("idmap_HGNC.txt", 2, 1, offset = 1).keys())

idmap_final = CIdMap()
idmap_final_ambiguous = CIdMap()



f = open("string_id_unique.txt")


for line in f.readlines():
    string_id = line.strip()
    ensp = string_id.split('.')[1]

    

f.close()




for ensp, hgncs in idmap_Biomart_Ensp_HGNC.iteritems():
    if len(hgncs) > 1:
        if idmap_Biomart_Ensp_Ensg.has_key(ensp):
            hgncs_from_HGNC = set()
            ensgs = idmap_Biomart_Ensp_Ensg.GetIds(ensp)
            for ensg in ensgs:
                if idmap_HGNC_Ensg_HGNC.has_key(ensg):
                    hgncs_from_HGNC.update(idmap_HGNC_Ensg_HGNC.GetIds(ensg))
                
            hgncs_both = hgncs.intersection(hgncs_from_HGNC)
            if len(hgncs_both) == 0:
                if manual_map.has_key(ensp):
                    idmap_final[ensp] = set([manual_map[ensp]])
                else:
                    if len(approved_symbol_names.intersection(hgncs)) == 1:
                        idmap_final[ensp] = approved_symbol_names.intersection(hgncs)
                    else:
                        print "No mapping in HGNC", ensp, ensgs, hgncs
                        idmap_final_ambiguous[ensp] = hgncs
            elif len(hgncs_both) > 1:
                if manual_map.has_key(ensp):
                    idmap_final[ensp] = set([manual_map[ensp]])
                else:
                    print ensp, hgncs_both
                    idmap_final_ambiguous[ensp] = hgncs_both
            else:
                idmap_final[ensp] = hgncs_both
                
        else:
            raise "No Ensg mapped!"
    else:
        idmap_final[ensp] = hgncs
        pass

fo = open("9606.alias_best.v9.0.17032014_2.tsv", 'w')
for ensp in idmap_final.keys():
    print >> fo, "9606\t%s\t%s\tHGNC symbol (Biomart, HGNC, Uniprot)" % (ensp, idmap_final.GetFirstId(ensp))


f = open("../data/9606.alias_best.v9.0.tsv")

for line in f.readlines():
    org, ensp, name, source = line.split()
    if not idmap_final.has_key(ensp):
        if source == "BioMart_HUGO" or source == "Ensembl_HGNC_curated_gene" or source == "Ensembl_HGNC":
            fo.write(line)
        elif idmap_final_ambiguous.has_key(ensp):
            print >> fo, "9606\t%s\t%s\tHGNC symbol (Biomart, HGNC, Uniprot)" % (ensp, idmap_final_ambiguous.GetFirstId(ensp))
        else:
            fo.write(line)

f.close()
fo.close()
