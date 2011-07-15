import MySQLdb

conn = MySQLdb.connect('ptilouis', 'biophenics', 'biophenics', 'biophenics')
cursor = conn.cursor()

platelist = None
platedir_to_plateid = {}

def digits(s):
    return filter(lambda c: c.isdigit(), s)

def get_plateid(platedir):
    global platelist
    if platelist is None:
        cursor.execute('SELECT PLATE_ID, PLATE_CODE FROM PLATE WHERE PLATE_CODE IS NOT NULL')
        plate_ids, plate_codes = zip(*cursor.fetchall())
        plate_codes = [digits(pc) for pc in plate_codes]
        platelist = [(pc, pid) for pc, pid in zip(plate_codes, plate_ids) if len(pc) >= 8]

    if platedir not in platedir_to_plateid:
        barcode = platedir
        # drop Incell postfix
        if barcode[-2] == '_' and barcode[-1] in '0123456789':
            barcode = barcode[:-2]
        # find last 8 digits
        barcode = digits(barcode)[-8:]
        for pc, pid in platelist:
            if barcode in pc:
                platedir_to_plateid[platedir] = pid
                break
        if platedir not in platedir_to_plateid:
            raise ValueError('Could not find plate from directory name %s' % (platedir))
    return platedir_to_plateid[platedir]

def lookup_genes(platedir, wellrow, wellcol):
    cursor.execute('''SELECT DUPLEX_NAME, DUPLEX_NUMBER FROM SI_RNA JOIN SI_RNA_WELL USING (SI_RNA_ID)
                      JOIN WELL USING (WELL_ID) JOIN PLATE USING (PLATE_ID) WHERE PLATE_ID=%d and WELL_ROW=%s and WELL_COL=%s''' %
                   (get_plateid(platedir), wellrow, wellcol))
    return ','.join(['%s#%d' % (g,n) for g, n in cursor.fetchall()])

def lookup_chemicals(platedir, wellrow, wellcol):
    cursor.execute('''SELECT CHEMICAL_NAME FROM CHEMICAL JOIN CHEMICAL_WELL USING (CHEMICAL_ID)
                      JOIN WELL USING (WELL_ID) JOIN PLATE USING (PLATE_ID) WHERE PLATE_ID=%d and WELL_ROW=%s and WELL_COL=%s''' %
                   (get_plateid(platedir), wellrow, wellcol))
    return ','.join(['%s' % (c[0]) for c in cursor.fetchall()])

def lookup_treatment(t):
    print t
    if t == 'Genes':
        return lookup_genes
    if t == 'Chemicals':
        return lookup_chemicals
    else:
        return None

if __name__ == '__main__':
    print lookup_genes('384_20X_Dapi_Alexa_cy3_2010-0445_1', 1, 1)
    print lookup_chemicals('384_20X_Dapi_Alexa_cy3_2010-0445_1', 1, 1)
    print get_plateid('384_20X_Dapi_Alexa_cy3_2010-0445_1')
