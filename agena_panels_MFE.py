import os
import json

'''
TM = {
    'min': 60,
    'max': 72
}

panel = {
    'AG57U': 'CACAGCTTGCTCACA',
    'AG59U': 'CCCATCATTGCAGACC',
    'AG60U': 'CAGGTGAGCTTTCAGT',
    'AG62U': 'AAAGGTAGCACACGAGG',
    'AG63U': 'TTGGCCTTTGACAAATAT'
}

data = {
    'new': {
        'seq': "GATCGTAGCTGATCGTAGTCGAT",
        'direction': "F",
        'snp': ["T", "C"],
        'weight': 5024
    },
    'new2': {
        'seq': "ATAGATAGGGGGATTGCTCACA",
        'direction': "R",
        'snp': ["T", "C"],
        'weight': 5024
    },
}
'''


def makeSeqUpperCase(panel):
    for key, value in panel.items():
        panel[key] = value.upper()

def makeMFEinput (data):
    MFEinput = ""
    for key, value in data.items():
        MFEinput+= ">{name}\n{seq}\n".format(name = key, seq = value.get("seq"))
    return MFEinput

def mfeInvoke (MFEinput):
    with open("temp/input.fa", "w") as file:
        file.write(MFEinput)
    bashCommand = "./mfeprimer dimer -i temp/input.fa --out temp/outputMFE --dg -3.5 --score 4 --mismatch 2 --mono 17 --diva 2 --dntp 0.2 --oligo 95 --json"
    os.system(bashCommand)
    path = "/home/alina/Desktop/GitHub/agena_panels/temp/outputMFE.json"
    while not os.path.exists(path):
        pass
    with open("temp/outputMFE.json") as file:
            readed = file.read()
    resultMFE = json.loads(readed)
    return resultMFE

def filterOne (resultMFE,numberOfPrimers,TM,data):
    dictCheckedUEP = dict()
    for i in range (0, numberOfPrimers):
        name = resultMFE['PrimerList'][i]['Seq']['ID']
        tm = resultMFE['PrimerList'][i]['Tm']
        dictCheckedUEP[name]=tm

    if resultMFE['DimerList'] != None:
        setOfDimers = set()
        for i in resultMFE['DimerList']:
            setOfDimers.add(i['S1']['ID'])
            setOfDimers.add(i['S2']['ID'])
    else: 
        setOfDimers = set("None")

    dictCheckedUEPafterFilter = dict()
    for key, value in dictCheckedUEP.items():
        if key not in setOfDimers and value > float(TM['min']) and value < float(TM['max']):
            data[key]["TM"] = "%.2f" % value #добавляем инфу про ТМ
            dictCheckedUEPafterFilter[key] = data[key]
    print("HAHA" + str(dictCheckedUEPafterFilter))
    return dictCheckedUEPafterFilter

def filterTwo (filter1, panel):
    acceptableUEP = dict()
    for key, value in filter1.items():
        panel[key] = value
        MFEinput = ""
        for key, value in panel.items():
            MFEinput+= ">{name}\n{seq}\n".format(name = key, seq = value)
        resultMFE2 = mfeInvoke(MFEinput)
        setOfDimers = set()
        if resultMFE2['DimerList'] != None:
            for i in resultMFE2['DimerList']:
                setOfDimers.add(i['S1']['ID'])
                setOfDimers.add(i['S2']['ID'])
            if key not in setOfDimers:
                acceptableUEP[key] = value
        else:
            acceptableUEP[key] = value
    return (acceptableUEP)

def addUEP(data, panel, TM):
    
    makeSeqUpperCase(panel)
    numberOfPrimers = len(data.keys())

    dir_temp = "/home/alina/Desktop/GitHub/agena_panels/temp"
    if not os.path.exists(dir_temp):
        os.makedirs(dir_temp)

    MFEinput = makeMFEinput (data)

    resultMFE1 = mfeInvoke(MFEinput)
    filter1 = filterOne(resultMFE1,numberOfPrimers,TM,data)
   # print("pass filter 1", filter1)

    filter2 = filterTwo(filter1, panel)
  #  print("pass filter 2", filter2)

    return filter2

#addUEP(data, panel, TM)

