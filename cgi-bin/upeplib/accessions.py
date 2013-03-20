import gzip

def seqtocaps(seq):
    temp = ''
    for i in seq:
        if ord(i) > 96:
            temp = temp + chr(ord(i)-32)
        else: temp
    return temp

def ACCloc(accession):
    filename = 'ACCcompletecompact.dict.rna.gbff'
    filepath = '../RefSeq/' + filename
    readfile = open(filepath,'rb')
    readfile.seek(0, 2)
    line = readfile.readline()
    low = 0
    high = int((readfile.tell()/80))
    pos = (low+high)/2
    while True:
        oldpos = pos
        readfile.seek(pos*80)
        line = readfile.readline()
        if accession < line.rsplit('\t')[0]:
            high = pos
            pos = (low+high)/2
        elif accession > line.rsplit('\t')[0]:
            low = pos
            pos = (low+high)/2
        else:
            return line.rstrip().split('\t')
        if pos == oldpos:
            return None
    readfile.close()

def GIloc(accession):
    compare = str(accession)
    filename = 'GIcompletecompact.dict.rna.gbff'
    filepath = '../RefSeq/' + filename
    readfile = open(filepath,'rb')
    readfile.seek(0, 2)
    line = readfile.readline()
    low = 0
    high = int((readfile.tell()/40))
    pos = (low+high)/2
    while True:
        oldpos = pos
        readfile.seek(pos*40)
        line = readfile.readline()
        compareto = line.rsplit('\t')[0]
        if compare < compareto:
            high = pos
            pos = (low+high)/2
        elif compare > compareto:
            low = pos
            pos = (low+high)/2
        else:
            return ACCloc(line.rsplit('\t')[1].rstrip())
        if pos == oldpos:
            return None
    readfile.close()    

def getCDS(accession):
    if len(accession) > 3:
        if accession[:3] == 'GI:':
            return GIloc(accession[3:])
        elif (accession[:3] == 'NM_' or accession[:3] == 'XM_'):
            return ACCloc(accession)
        else:
            return None
    else:
        return None

def getmRNA(accession, minsize, maxsize, gracelength):
    takejoins = True
    details = getCDS(accession)
    if details:
        handle = gzip.open('../RefSeq/' + details[2],'rb','9')
        try:
            handle.seek(int(details[1]))
            while True:
                line = handle.readline()
                if line[0:10] == 'DEFINITION':
                    definition = line[12:].rstrip()
                    line = handle.readline()
                    while True:
                        if line[0] == ' ':
                            definition = definition + ' ' + line[12:].rstrip()
                            line = handle.readline()
                        else:
                            break                        
                if line[5:8] == 'CDS':   
                    if line[21:25] == 'join':
                        if takejoins:
                            more = returnjoins(line[25:].rstrip())
                            CDSstart = int(more[0][0])
                            CDSend = int(more[len(more)-1][1])
                    else:
                        i = line[21:]
                        for j in range(0,len(i)):
                            if i[j] == '.':
                                CDSstart = int(i[0:j])
                                CDSend = int(i[j+2:].rstrip())
                                break                 
                elif line[0:6] == 'ORIGIN':
                    temp = ''
                    while True:
                        line = handle.readline()
                        if not line or line[0:2] == '//':
                            break
                        for i in line[10:].rstrip():
                            if not i == ' ':
                                temp = temp + i                         
                    uORFs = []
                    if maxsize:
                        for i in range(0, CDSstart):
                            if temp[i] == 'a' and temp[i+1] == 't' and temp[i+2] == 'g':
                                j = 0
                                while (i + j) < CDSstart + gracelength:
                                    j = j + 3
                                    if temp[i+j] == 't' and (temp[i+j+1] == 'a' and (temp[i+j+2] == 'g' \
                                        or temp[i+j+2] == 'a') or (temp[i+j+1] == 'g' and (temp[i+j+2] == 'a'))):
                                        uORFs.append([i,i+j])
                                        break    
                    if uORFs:
                        j = 0
                        while True:
                            if uORFs[j][1]-uORFs[j][0] < minsize or uORFs[j][1]-uORFs[j][0] > maxsize:
                                del uORFs[j]
                                j = j - 1
                            j = j + 1
                            if j == len(uORFs):
                                break
                    temp = seqtocaps(temp)               
                    return [temp, [CDSstart, CDSend], uORFs, definition]       
        finally:
            handle.close()                           
    return None





