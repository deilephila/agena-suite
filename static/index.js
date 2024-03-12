const inputFile = document.getElementById("inputFile");
function handleFile(file) {
  //console.log("file:", file);
    let reader = new FileReader();
    reader.onload = () => {finalOutput = regexSearch(reader.result)};
    reader.readAsText(file);
}

function regexSearch(text) {
        const regex = /W\d[^W]*/g
        const matched = [... text.matchAll(regex)];
        fromFile = matched;
        let arrayOfObjects = renderResult(matched);
        return arrayOfObjects
    }

function makeAllele (oneRS) {
    let array = [];
    for (let i = 14; i<25; i += 3) {
        let ntpToAdd = oneRS[i+2].split("").at(-1);
        let snp = oneRS[11] == "F" ? ntpToAdd : cN[ntpToAdd];
        let obj = {
            allele: oneRS[i],
            EPmass: oneRS[i+1],
            EPseq: oneRS[i+2]
        }
        array.push(obj)
    }
    return array
}

function renderResult(res){
    let arrOfObjFromFile = [];
    for (let i=0; i<res.length; i++) {
        let oneRS = res[i][0];
        oneRS = oneRS.split("\t");
        let objFromFile = {
            name: oneRS[2],
            dir: oneRS[11],
            UEPmass: oneRS[12],
            seq: oneRS[13],
            EP: makeAllele(oneRS)
        }
        console.log(objFromFile)
        let snp = [];
        objFromFile.EP.forEach((e)=> {
            e.allele && snp.push(e.allele)
        }); 
        objFromFile.snp = snp;    
        arrOfObjFromFile.push(objFromFile);
    }
    console.log(arrOfObjFromFile);
    let arrayOfObjects = fromXLStoInput (arrOfObjFromFile);
    return arrayOfObjects;
}

function fromXLStoInput (arrOfObjFromFile) {
    let arrayOfObjects = [];
    arrOfObjFromFile.forEach((rs) => {
        let objectOfUEP = {
            name: rs.name,
            seq: rs.seq.toUpperCase(),
            dir: rs.dir,
            snp: rs.snp,
            UEPmass: calculateUEPMass(rs.seq.toUpperCase()),
        }
        objectOfUEP.EP = calculateEPMass(objectOfUEP);
        arrayOfObjects.push(objectOfUEP);
    })
    console.log(arrayOfObjects);
    
    let arrayConflicts = massCompare(arrayOfObjects);
    console.log(arrayConflicts);
    let htmlStatus = CheckSetOutputStatus(arrayConflicts);
    document.getElementById("outputStatus").innerHTML = htmlStatus;
    let htmlTable = CheckSetOutputTable (arrayOfObjects);
    document.getElementById("outputSetOfUEP").innerHTML = htmlTable; //рендер таблицы с данными о панели
    makeVisible("invisDownload") //после рендера результатов должна появиться кнопка Download
    return arrayOfObjects
}

//базовая функция - для определения выбранного radio button - универсальная! на вход - name radio button ("r1")
function radioButton(name) {
    let rad=document.getElementsByName(name);
    let output;
    for (let i=0; i<rad.length; i++) {
        if (rad[i].checked) {
            output = i;
        }
    }
    return output
}

// изменение value внутри инпута для массы каждой буквы в зависимости от выбранного ТМ (radio button)
function chooseTM() {
    let outputFromR2 = radioButton("r2");
    let TMiPlex = {
        'A': 271.2,
        'T': 327.1,
        'G': 287.3,
        'C': 247.2
    }
    let TMddNTP = {
        'A': 296.2,
        'T': 287.1,
        'G': 312.3,
        'C': 272.2
    }
    let TM3M = {
        'A': 601.4,
        'T': 287.1,
        'G': 618,
        'C': 578.5
    }
    let TM1M = {
        'A': 296.2,
        'T': 577.8,
        'G': 312.3,
        'C': 272.2
    }
    let TMU = {
        'A': 296.2,
        'T': 441,
        'G': 312.3,
        'C': 272.2
    }


    let arrayOfTM = [TMiPlex, TMddNTP, TM3M, TM1M, TMU];
    document.getElementById("InputMassA").value = arrayOfTM[outputFromR2]['A'];
    document.getElementById("InputMassT").value = arrayOfTM[outputFromR2]['T'];
    document.getElementById("InputMassC").value = arrayOfTM[outputFromR2]['C'];
    document.getElementById("InputMassG").value = arrayOfTM[outputFromR2]['G'];
    document.getElementById("InputCheckMass").dispatchEvent(new Event("input"));
    document.getElementById("InputSetOfUEP").dispatchEvent(new Event("input"));
}

// функция для записи значений масс терминаторов (введеных пользователем), output - array
function inputTerminationMass () {
    let termMasses = {};
    termMasses["A"] = parseFloat(document.getElementById('InputMassA').value.replace(',', '.'));
    termMasses["T"] = parseFloat(document.getElementById('InputMassT').value.replace(',', '.'));
    termMasses["C"] = parseFloat(document.getElementById('InputMassC').value.replace(',', '.'));
    termMasses["G"] = parseFloat(document.getElementById('InputMassG').value.replace(',', '.'));
    return (termMasses);
}

//базовая функция - расчет массы UEP
function calculateUEPMass(UEPsequence) {
    let A = 0; let G = 0; let C = 0; let T = 0;
    let seq = UEPsequence;
    let UEPLength = seq.length;
    for (i in seq) {
        if (seq[i] == 'A') {A += 1;}
        else if (seq[i] == 'G') {G += 1;}
        else if (seq[i] == 'C') {C += 1;}
        else if (seq[i] == 'T') {T += 1;}
        else {throw new Error("invalid character in sequence")}}
    let UEPmass = Number((249.2*A + 225.2*C + 265.2*G + 240.2*T + 64*(UEPLength - 1) + 2.95).toFixed(1));
    return (UEPmass)};

//мапка для комплементарных букв    
let cN = {
    'C': 'G',
    'A': 'T',
    'G': 'C',
    'T': 'A',
}

//базовая функция расчета массы EP
function calculateEPMass (objectOfUEP) {
    let terminationMasses = inputTerminationMass();
    let mass = objectOfUEP.UEPmass;
    let snp = objectOfUEP.snp; // на вход - ["C", "T"] - надо!
    let dir = objectOfUEP.dir; // на вход - 'F', формат string
    let seq = objectOfUEP.seq;
    let arrayOfObjForEP = [];
    
    for (i in snp) {
        let termNTP = dir == 'F' ? snp[i] : cN[snp[i]];

        let EPobj = {
            allele: snp[i],
            EPmass: Math.round((mass + terminationMasses[termNTP])*10)/10,
            EPseq: seq + termNTP
        };
        arrayOfObjForEP[i] = EPobj;
     }
    return (arrayOfObjForEP)
 }

//checkSet func 0 - проверка панели на наличие конфликтов, начинает работать при инпуте в textarea
function checkSet(text) {
    if (!text) return;
    let outputFromR1 = radioButton("r1"); //базовая функция для определения формата ввода данных - seq или mass (сейчас по умолчанию всегда seq!)
    let strings = text.trim().split("\n").toString().split(","); //чтение инпута
    let arrayOfObjects = []; //создаем будущий аутпут (инфа о всей панели)
    
    strings.forEach(str => { // каждая строка - один UEP, создаем для него obj
        str = str.split("	");
        let objectOfUEP = {
            name: str[0],
            seq: outputFromR1 == 0
                ? str[1].toString()
                : undefined,
            dir: str[2],
            snp: str[3].split("/"), //убираем слэш, нам нужен формат [T,C]
            UEPmass: outputFromR1 == 0
                ? calculateUEPMass(str[1].toString().toUpperCase()) //считаем массу, если на вход - seq
                : str[1],
        }
        objectOfUEP.EP = calculateEPMass(objectOfUEP); //добавляем инфу об EP с помощью базовой функции
        arrayOfObjects.push(objectOfUEP); //пушим объект по одному UEP в общий array (будущий output)
    })
    let arrayConflicts = massCompare(arrayOfObjects); //базовая функция по проверке конфликтов
    let htmlStatus = CheckSetOutputStatus(arrayConflicts); // checkSet func 1 - формирование статуса конфлитов (есть или нет)
    document.getElementById("outputStatus").innerHTML = htmlStatus; // рендер статуса конфликтов 
    let htmlTable = CheckSetOutputTable (arrayOfObjects); // checkSet func 2 - формирование таблицы c данными о панели
    document.getElementById("outputSetOfUEP").innerHTML = htmlTable; //рендер таблицы с данными о панели
    makeVisible("invisDownload") //после рендера результатов должна появиться кнопка Download
    finalOutput =  arrayOfObjects;
    return arrayOfObjects
}

//базовая функция для сравнения масс UEP и EP
function massCompare(arrayOfObjects) {
    let arrayOfMasses = [];
    arrayOfObjects.forEach((oneUEP) => {
        let nameUEP = oneUEP.name
        let objUEP = {
            name: nameUEP,
            mass: oneUEP.UEPmass
        };
        arrayOfMasses.push(objUEP);
        oneUEP.EP.forEach((allele) => {
            let objEP = {
            name: nameUEP,
            mass: allele.EPmass
        }
        arrayOfMasses.push(objEP);
        })
    })
    let arrayConflicts = [];
    for (let i=0; i < arrayOfMasses.length; i++) {
        let mass1 = arrayOfMasses[i].mass;
        let name1 = arrayOfMasses[i].name;
        for (let k = i + 1; k < arrayOfMasses.length; k++) {
            if (arrayOfMasses[k].name != name1) {
                let mass2 = arrayOfMasses[k].mass;
                let delta = mass2 - mass1;
                if (delta < 30 && delta > -30) {
                    arrayConflicts.push([arrayOfMasses[i], arrayOfMasses[k]]);
                }}}}
    let names = [];
    arrayConflicts.forEach((array) => {
        let str = [array[0].name, "/*!/", array[1].name].join("");
        names.push(str);
    });
    uniqueArray = names.filter(function(item, pos) {
        return names.indexOf(item) == pos;
    })
    let uniqueArrayFinal = [];
    uniqueArray.forEach((conflict) => {
        uniqueArrayFinal.push(conflict.split("/*!/"))
    })
    return (uniqueArrayFinal)
}

//checkSet func 1 - формирование аутпута для статуса проверки панели (есть конфликты или нет)
function CheckSetOutputStatus (uniqueArrayFinal) {
    let html = "";
    if (uniqueArrayFinal.length == 0) {
        html = "Great UEPs!";
        makeVisible("invis2")
    } else {
        uniqueArrayFinal.forEach((conflict) => {
            html += `<div> Attention! There is a conflict between ${conflict[0]} and ${conflict[1]}! </div>`;
        })
    }
    return html
}

//checkSet func 2 - формирование аутпута таблицы с данными по каждому олигу
function CheckSetOutputTable (arrayOfObjects) {
    let html = `<tr>
        <th>Name</th>
        <th>Sequence</th>
        <th>Direction</th>
        <th>SNP</th>
        <th>UEP mass</th>
        </tr>`
    arrayOfObjects.forEach((array) => {
        let snp = array.snp.join("/");
        html += `<tr>
                    <td>${array.name}</td>
                    <td>${array.seq}</td>
                    <td>${array.dir}</td>
                    <td>${snp}</td>
                    <td>${array.UEPmass}</td>
                </tr>`;
    })
    return html
}

function CheckMassOutputTable (arrayOfObjects) {
    let html = `<tr>
        <th>Sequence</th>
        <th>UEP mass</th>
        <th>+A</th>
        <th>+T</th>
        <th>+C</th>
        <th>+G</th>
        </tr>`
    arrayOfObjects.forEach((array) => {
        html += `<tr>
            <td>${array.seq}</td>
            <td>${array.UEPmass}</td>
            <td>${array.EP[0].EPmass}</td>
            <td>${array.EP[1].EPmass}</td>
            <td>${array.EP[2].EPmass}</td>
            <td>${array.EP[3].EPmass}</td>
            </tr>`;
    })
    document.getElementById("outputCheckMass").innerHTML = html;
}

function checkMass (text) {
    if (!text) return;
    let strings = text.trim().split("\n").toString().split(",");
    let arrayOfObjects = [];
    strings.forEach(str => {
        let objectOfUEP = {
            seq: str.toString(),
            UEPmass: calculateUEPMass(str.toString().toUpperCase()),
            dir: "F",
            snp: ["A", "T", "C", "G"]
        }
        objectOfUEP.EP = calculateEPMass(objectOfUEP);
        arrayOfObjects.push(objectOfUEP);
    })
    console.log(arrayOfObjects);
    CheckMassOutputTable(arrayOfObjects);
}

document.getElementById("InputCheckMass").addEventListener("input", () => checkMass(document.getElementById("InputCheckMass").value));
document.getElementById("InputSetOfUEP").addEventListener("input", () => checkSet(document.getElementById("InputSetOfUEP").value));
document.getElementsByName('r1').forEach( r =>
    r.addEventListener("input", () => checkSet(document.getElementById("InputSetOfUEP").value))
)
inputFile.addEventListener("change",(e)=>handleFile(e.target.files[0]));

/* при переключении ТМ (radio button) изменяются value внутри инпутов масс букв*/
document.getElementsByName('r2').forEach( r =>
    r.addEventListener("change", () => chooseTM()));

document.getElementById("buttonDownload").addEventListener("click", Download)

function ConvertToCSV(objArray) {
  let array = typeof objArray != 'object' ? JSON.parse(objArray) : objArray;
  let str = '';
  for (let i = 0; i < array.length; i++) {
    var line = '';
    for (var index in array[i]) {
      if (line != '') {
        line += ','
        }
      line += array[i][index];
    }
    str += line + '\r\n';
  }
  return str;
}

function makeJSON (finalOutput, fromFile) {
    let json = [];
    let firstString = {
        "WELL": "WELL", "TERM": "TERM", "SNP_ID": "SNP_ID",
        "2nd-PCRP": "2nd-PCRP", "1st-PCRP": "1st-PCRP", "AMP_LEN": "AMP_LEN",
        "UP_CONF": "UP_CONF", "MP_CONF": "MP_CONF", "Tm(NN)": "Tm(NN)",
        "PcGC": "PcGC", "PWARN": "PWARN", "UEP_DIR": "UEP_DIR", "UEP_MASS": "UEP_MASS", 
        "UEP_SEQ": "UEP_SEQ", "EXT1_CALL": "EXT1_CALL", "EXT1_MASS": "EXT1_MASS",
        "EXT1_SEQ": "EXT1_SEQ", "EXT2_CALL": "EXT2_CALL", "EXT2_MASS": "EXT2_MASS",
        "EXT2_SEQ": "EXT2_SEQ", "EXT3_CALL": "EXT3_CALL", "EXT3_MASS": "EXT3_MASS",
        "EXT3_SEQ": "EXT3_SEQ", "EXT4_CALL": "EXT4_CALL", "EXT4_MASS": "EXT4_MASS",
        "EXT4_SEQ": "EXT4_SEQ", "1stPAUSE": "1stPAUSE"
    };
    json.push(firstString);
    if (fromFile[0] == undefined) {
        console.log(fromFile[0] == undefined)
        finalOutput.forEach((oneRS) => {
            let string = {
            "WELL": " ", "TERM": "iPlex", "SNP_ID": oneRS.name,
            "2nd-PCRP": " ", "1st-PCRP": " ", "AMP_LEN": " ",
            "UP_CONF": " ", "MP_CONF": " ", "Tm(NN)": " ",
            "PcGC": " ", "PWARN": " ", "UEP_DIR": oneRS.dir, "UEP_MASS": oneRS.UEPmass, 
            "UEP_SEQ": oneRS.seq,
            "EXT1_CALL": oneRS.EP[0].allele, "EXT1_MASS": oneRS.EP[0].EPmass, "EXT1_SEQ": oneRS.EP[0].EPseq,
            "EXT2_CALL": oneRS.EP[1] == undefined ? "" : oneRS.EP[1].allele, "EXT2_MASS": oneRS.EP[1] == undefined ? "" : oneRS.EP[1].EPmass, "EXT2_SEQ": oneRS.EP[1] == undefined ? "" : oneRS.EP[1].EPseq,
            "EXT3_CALL": oneRS.EP[2] == undefined ? "" : oneRS.EP[2].allele, "EXT3_MASS": oneRS.EP[2] == undefined ? "" : oneRS.EP[2].EPmass, "EXT3_SEQ": oneRS.EP[2] == undefined ? "" : oneRS.EP[2].EPseq,
            "EXT4_CALL": oneRS.EP[3] == undefined ? "" : oneRS.EP[3].allele, "EXT4_MASS": oneRS.EP[3] == undefined ? "" : oneRS.EP[3].EPmass, "EXT4_SEQ": oneRS.EP[3] == undefined ? "" : oneRS.EP[3].EPseq,
            "1stPAUSE": " "};
            json.push(string);
        })
        console.log(json);
    } else {
        for (let i=0; i<fromFile.length; i++) {
            let oneRS = fromFile[i][0];
            oneRS = oneRS.split("\t");
            console.log(oneRS);
            let string = {
                "WELL": oneRS[0], "TERM": oneRS[1], "SNP_ID": oneRS[2],
                "2nd-PCRP": oneRS[3], "1st-PCRP": oneRS[4], "AMP_LEN": oneRS[5],
                "UP_CONF": oneRS[6], "MP_CONF": oneRS[7], "Tm(NN)": oneRS[8],
                "PcGC": oneRS[9], "PWARN": oneRS[10], "UEP_DIR": oneRS[11], "UEP_MASS": finalOutput[i].UEPmass, 
                "UEP_SEQ": oneRS[13], "EXT1_CALL": oneRS[14], "EXT1_MASS": finalOutput[i].EP[0].EPmass,
                "EXT1_SEQ": oneRS[16], "EXT2_CALL": oneRS[17], "EXT2_MASS": finalOutput[i].EP[1].EPmass,
                "EXT2_SEQ": oneRS[19], "EXT3_CALL": oneRS[20], "EXT3_MASS": finalOutput[i].EP[2] == undefined ? "" : finalOutput[i].EP[2].EPmass,
                "EXT3_SEQ": oneRS[22], "EXT4_CALL": oneRS[23], "EXT4_MASS": finalOutput[i].EP[3] == undefined ? "" : finalOutput[i].EP[3].EPmass,
                "EXT4_SEQ": oneRS[25], "1stPAUSE": " ",}
            json.push(string);
        }
    }
    return json
}

function Download() {
    console.log(fromFile)
    let json = fromFile[0] == undefined ? finalOutput : fromFile;
    let newjson = makeJSON (finalOutput, fromFile);
    let csv = ConvertToCSV(newjson);
    let downloadLink = document.createElement("a");
    let blob = new Blob(["\ufeff", csv]);
    let url = URL.createObjectURL(blob);
    downloadLink.href = url;
    downloadLink.download = "data.csv";

    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
};

//далее модуль AddUEP

document.getElementById("buttonAddUEP").addEventListener("click", () => makeVisible("invis")); //нажатие AddUEP -> инпуты для этого модуля становятся видимыми (базовая ф-ия makeVisible)

//класс WC - для рендера аутпутной таблицы модуля AddUEP
class WCoutTable extends HTMLElement {
    shadowDOM = this.attachShadow({mode: "closed"});
    constructor() {
        super();
        const template = document.getElementById("template1").content.cloneNode(true);
        this.shadowDOM.appendChild(template);
    }
    set html(htmlOutput) {
        this.shadowDOM.querySelector("table").innerHTML = htmlOutput
    }
    set warning(htmlWarning) {
        this.shadowDOM.querySelector("div").innerHTML = htmlWarning
    }
}

function clear(elementName) {
        while (elementName.firstChild) {
            elementName.firstChild.remove();
        }
    }

const out = document.getElementById("out");
customElements.define("table-out", WCoutTable);

const outAddUEP = new WCoutTable(); //создаем константу, которую потом будем изменять (в нее будет идти разный аутпут)
document.body.appendChild(outAddUEP); //добавляем эту константу в ДОМ


//базовая функция makeVisible - делает элементы класса (имя класса подается на вход) видимыми

function makeVisible(className) {
    document.getElementsByClassName(className)[0].style= "visibility: visible";
}

document.getElementById("addOligo").addEventListener("click", () => addUEP()); //нажатие кнопки Add (после ввода данных о added UEP) -> функция AddUEP (модуль)

// AddUEP func 0 - главная функция по добавлению UEP в существующий assay, вызывается кнопкой, результат - вывод таблицы
function addUEP () {
    let ampseq = document.getElementById("amp_oligoToAdd").value; //получение инпутных данных о посл-ти ампликона
 //   let dir = radioButton("r3")==0 ? "F" : "R"; //получение инпутных данных о направлении UEP, базовая ф-ия
 //   ampseq = dir=="F" ? ampseq : reverseComplement(ampseq); // AddUEP func 4 - если dir = R, то amp -> reverse complement
    let matched = readAmpSeq (ampseq); //AddUEP func 1 - прочтение посл-ти, разделение на forw, rev, snp
    let forward = matched[0].toString(); //посл-ть левой части ампликона (для forward UEP)
    let reverse = matched[2].toString(); //посл-ть правой части ампликона (для reverse UEP)
    let snp = matched[1].toString().split("/"); //определяем SNP, сразу переводим в формат ["A", "C"]
    let forUEPmax = cutAmp(forward); //AddUEP func 3 - cutAmp - обрезает ампликон forward до размера <9000 Да
    let revUEPmax = cutAmp(reverseComplement(reverse)); //AddUEP func 2 и 3 - аналогично, но т.к. rev, сначала надо сделать reverse complement (func 2)
    let forSuitableUEP = chooseUEP(forUEPmax, "F", snp); //AddUEP func 4 - chooseUEP - аутпут - array из подходящих UEP (не вступают в конфликт с сущ. assay) 
    let revSuitableUEP = chooseUEP(revUEPmax, "R", snp); //AddUEP func 4 - chooseUEP - аналогично для reverse
    let suitableUEP = forSuitableUEP; //создаем общий массив suitable UEP, в который сразу включаем forward UEP
    revSuitableUEP.forEach((i) => { //а затем с помощью этой функции пушим туда же все reverse UEP
        suitableUEP.push(i)
    })
    //если suitableUEP пустой, то не надо даже делать request, надо сразу выводить warning
    if (suitableUEP.length == 0) {
        let html = "Невозможно подобрать UEP подходящей массы :(";
      //  clear(out);
        outAddUEP.warning = html 
        return
    }

    //тут создаем набор UEP на проверку на бэке
    let suitableUEP_forBack = new Object
    suitableUEP.forEach ((i) => {
        suitableUEP_forBack[i.name] = {
            "seq": i.seq,
            "direction": i.dir,
            "snp": i.snp,
            "weight": i.UEPmass
        } 
    })

    //далее надо создать панель в нужном для бэка формате
    let arrayOfObjects = checkSet(document.getElementById("InputSetOfUEP").value);
    let panelBack = new Object;
    arrayOfObjects.forEach((i) => {
        panelBack[i.name] = i.seq
    })

    //тут создаем значения TM в нужном для бэка формате
    let TMBackend = {
        'min': 60,
        'max': 72
    }

    //а тут формируем собственно object для передачи на бэк
    let forBackend = {
        "UEP": suitableUEP_forBack,
        "Panel": panelBack,
        "TM": TMBackend
    }
    console.log(forBackend)
    suitableUEP_request (forBackend)
}


function suitableUEP_request (forBackend) {
    async function postAddUEP() {
        const response = await fetch('http://127.0.0.1:8000/addUEP', {
            method: "POST",
            headers: {
                "Content-Type": "application/json; charset=utf-8"
            },
            body: JSON.stringify(forBackend)
        });
        console.log(await response)
        const result = await response.json();
        console.log("result in async function", result)

        let htmlOtput = addUEPoutputNEW (result);
      //  clear(out); 
        outAddUEP.html = htmlOtput; 
        return result
    }
    postAddUEP()
}



//Add UEP func 1 - чтение инпута регексом,на выход - array из трех мэтчей (forward, reverse, snp)
function readAmpSeq (ampseq) {
    const regex = /[A-Z\/]{1,}/gm; 
    const matched = [... ampseq.matchAll(regex)];
    return matched;
}

// AddUEP func 2 - функция для превращения reverse amp в forward amp
function reverseComplement (ampseq) {
    let complementAmpseq = "";
    //с помощью мапки превращаем буквы в комплеметарные
    for (let i=0; i<ampseq.length; i++) {
        complementAmpseq += cN[ampseq[i]]
    }
    let reverseComplementAmpseq = complementAmpseq.split("").reverse().join("") //делаем реверс цепи
    return reverseComplementAmpseq
}

// AddUEP func 3 - функция для обрезки ампликона до размера <9000 Да
function cutAmp(ampseq) {
    let maxMass = 0; //масса олига (конечная) на каждом цикле (+64*(length...))
    let UEPmax = ""; //последовательность олига, растет на каждом цикле
    let UEPmassCount = 0; //масса олига неконечная (!), масса только букв (без +64...), нужна для счета
    for (let i=ampseq.length-1; i>=0; i--) { //здесь - forward UEP, т.е. последняя буква (-1) - соседняя с SNP, поэтому надо удалять нуклеотиды с начала (0) до тех пор, пока масса не будет < 9000
        if (maxMass > 9000) break
        if (ampseq[i] == 'A') {
            maxMass = Number((Number(UEPmassCount)+249.2 + 64*(UEPmax.length - 1) + 2.95).toFixed(1));
            UEPmassCount += 249.2;
            UEPmax += "A";
        } else if (ampseq[i] == 'G') {
            maxMass = Number((Number(UEPmassCount)+265.2 + 64*(UEPmax.length - 1) + 2.95).toFixed(1));
            UEPmassCount += 265.2;
            UEPmax += "G";
        } else if (ampseq[i] == 'C') {
            maxMass =Number((Number(UEPmassCount)+225.2 + 64*(UEPmax.length - 1) + 2.95).toFixed(1));
            UEPmassCount += 225.2;
            UEPmax += "C";
        } else if (ampseq[i] == 'T') {
            maxMass = Number((Number(UEPmassCount)+240.2 + 64*(UEPmax.length - 1) + 2.95).toFixed(1));
            UEPmassCount += 240.2;
            UEPmax += "T";
        } else {throw new Error("invalid character in sequence")}
    }
    //мы получили reverse amp, надо его перевернуть и убрать первую букву
    UEPmax = UEPmax.slice(0,-1).split("").reverse().join("");
    return UEPmax
};

// addUEP func 4 - подбор подходящих вариантов UEP, запуск из func 0
function chooseUEP(UEPmax, dir, snp) {
    let suitableUEP = []; //создаем пустой массив, который потом будет аутпутом
    let UEPseqVar = UEPmax; //временная последовательность UEP, которая будет уменьшаться на 1 нуклеотид в каждом цикле
    let UEPmassVar = calculateUEPMass(UEPseqVar); //расчет массы временной посл-ти UEP через базовую ф-ию 
    //цикл, в котором на каждой итерации создается obj для UEP, проверка конфликтов и расчет массы (>4300 и <9000)
    let i = 0; //создаем счетчик (только ради нумерации в названии new UEP)
    while (UEPmassVar > 4300) {
        //создаем для каждого UEP object
        let newUEP = {
            name: `new${i}`,
            seq: UEPseqVar,
            dir: dir,
            snp: snp,
            UEPmass: UEPmassVar
        }
        newUEP.EP = calculateEPMass(newUEP); //добавляем в только что созданный obj инфу об EP (базовая ф-ия)
        //берем актуальный arrayOfObjects (текущий набор праймеров, загруженный пользователем)
        let arrayOfObjects = checkSet(document.getElementById("InputSetOfUEP").value);
        //исключаем UEP, у которых масса EP > 9000 (добавляем в arrayOfObject только те UEP, EP которых < 9000)
        let lessThan9000 = 0;
        newUEP.EP.forEach((i)=> {
            if (i.EPmass > 9000) {
                lessThan9000 += 1;
            }
        })
        if (lessThan9000 == 0) {
            arrayOfObjects.push(newUEP);
            //и тут же проверяем обновленный arrayOfObjects на конфликты, если их нет, то новый UEP пушится в suitableUEP
            if (massCompare(arrayOfObjects).length == 0) {
                suitableUEP.push(newUEP); //этот array - будущий аутпут
            }
        }
        UEPseqVar = UEPseqVar.slice(1); //убираем первую букву UEPseqVar, чтобы в следующем цикле прогонять новую посл-ть 
        UEPmassVar = calculateUEPMass(UEPseqVar) //и сразу считаем массу нового UEPseqVar
        i ++ //и не забываем увеличить счетчик
    }
    return suitableUEP
}

// вывод результата подбора UEP в виде таблички - после MFE!
function addUEPoutputNEW (suitableUEP) {
    if (Object.keys(suitableUEP).length == 0) {
        let html = "Ни один UEP не удовлетворяет условиям мультиплекса :(";
      //  clear(out);
        outAddUEP.warning = html
        return
    }
    let html = `<tr>
        <th>Name</th>
        <th>Sequence</th>
        <th>Direction</th>
        <th>SNP</th>
        <th>UEP mass</th>
        <th>Tm</th>
        </tr>`
    let i = 1; //ставим тут счетчик, чтобы можно было включить порядковый номер праймера в название
    
    Object.keys(suitableUEP).forEach( (key) => {
        let snp = suitableUEP[key].snp.join("/");

        html += `<tr>
                <td>${key}-${suitableUEP[key].direction}</td>
                <td>${suitableUEP[key].seq}</td>
                <td>${suitableUEP[key].direction}</td>
                <td>${snp}</td>
                <td>${suitableUEP[key].weight}</td>
                <td>${suitableUEP[key].TM}</td>
            </tr>`;
        i++;
    });
    return html
}
