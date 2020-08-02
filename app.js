// Main function to calculate rho app

function firtsSubscript(f, thicks, rhos, dists, coilsType){
    const minRho = Math.min(...rhos)
    if (coilsType === 'H'){
        const logThick = Math.log(1000*thicks[0]/dists)  
        let KS = Math.floor(4*logThick)
        if (KS > 22){
            KS = 22
        } 
        else if(KS < 0){
            KS = 0
        }

        return KS

    }
    const logThick = Math.log(1E9*minRho/(f*dists**2))
    let KS = Math.floor(1.3*logThick)
    if (KS > 11){
        KS = 11
    }
    else if(KS < 0){
        KS = 0
    }
    return KS
}

function abcissa(KS, dists, coilsType){
    let K = 38
    let No = 6.10113727
    if (coilsType === 'H'){
        K = 51
        No = 8.75198087
    }

    const Y = new Array(K).fill(0)
    for (k=KS; k < K; k++){
        const temp1 = No - (k)*0.1*Math.log(10)
        const temp2 = Math.log(dists)
        Y[k] = (temp2 - temp1)
    }
    return [Y, K]
}


function Ro(rhos, f, lamda, thicks){
    const MU = (4E-7)*Math.PI

    const layer = rhos.length
    const V = new Array(layer).fill(0)
    for (j4=0; j4<layer; j4++){
        const gamma2 = 2*Math.PI*MU*f/rhos[j4]
        const tempc1 = math.complex(lamda ** 2, gamma2)

        V[j4] = math.sqrt(tempc1)
        
    }

    let P = 0
    let L = 0
    let M = 0
    let R = 0

    let I = rhos.length
    while (I > 0){
        I = I - 1
        if (I === 0){
            P = div(math.subtract( lamda, V[0]), math.add(lamda, V[0])) 
        }
        else{
            P = div(math.subtract( V[I-1], V[I]), math.add(V[I-1], V[I]))
            
        }
        
        L = math.complex(-2*thicks[I]*V[I].re, -2*thicks[I]*V[I].im)
        
        if (L.re > -50){
            M = math.exp(L)
            
        }
        else{
            M = 0
        }
        R = div(math.add(P, math.multiply(R, M)), math.add(1, math.multiply(P, R, M)))
              
    }
    
    return [R, P]
}

function mutualCouplingHomogen(dists, f, rhos, coilsType){
    const MU = (4E-7)*Math.PI
    
    const tempc1 = math.sqrt(math.complex(0,(dists ** 2) * MU * 2 * Math.PI*f / rhos[0]))

    const tempc2 = math.add(9, math.multiply(9, tempc1), math.multiply(4, math.pow(tempc1, 2)), math.pow(tempc1, 3) )
    
    const tempc3 =  math.multiply(tempc2, math.exp(math.multiply(-1,tempc1)))

    let Z = 0
    if (coilsType === 'H'){
        Z = div(math.multiply(2, math.subtract(9, tempc3)), math.pow(tempc1, 2))
    }
    return Z
}


function condApp(f, thicks, rhos, dists, coilsType){
    const MU = (4E-7)*Math.PI

    const ck = {
            'H' : [
                -0.00001787,  0.00000935,  0.00000375, -0.00001754,  0.00001084,
                -0.00000942,  0.00000456,  0.00000394, -0.00001576,  0.00003025,
                -0.00004683,  0.00006539, -0.00008669,  0.00011278, -0.00014748,
                 0.00019692, -0.00027055,  0.00038337, -0.00056557,  0.00085297,
                -0.00134318,  0.00224120, -0.00404751,  0.00812962, -0.01859531,
                 0.04821827, -0.13070863,  0.31328618, -0.51302191,  0.31003396,
                 0.34216522, -0.20142842, -0.36288158, -0.22914055, -0.03202792,
                 0.10252302,  0.16941035,  0.18559086,  0.17656063,  0.15523408,
                 0.13149777,  0.10841834,  0.08826593,  0.07109834,  0.05708674,
                 0.04555925,  0.03635105,  0.02892500,  0.02289634,  0.01843169,
                 0.01454755],
            'V' : [
                -0.00001323, 0.00003397, -0.00006292,  0.00010397,
                -0.00016337, 0.00025153, -0.00038653,  0.00060158,
                -0.00096163, 0.00160407, -0.00284834,  0.00552295,
                -0.01202839, 0.02983246, -0.08103193,  0.21267387,
                -0.43674023, 0.49063145,  0.04195061, -0.43651651,
                -0.22404767, 0.09474711,  0.26322713,  0.28286168,
                 0.23816634, 0.17652451,  0.12333557,  0.08243897,
                 0.05416555, 0.03489808,  0.02239027,  0.01424526,
                 0.00903156, 0.00574492,  0.00359571,  0.00231841,
                 0.00141123, 0.00094868
                ]
    }

    const KS = firtsSubscript(f, thicks, rhos, dists, coilsType)
    const [Y, K] = abcissa(KS, dists, coilsType)
    let Z =  mutualCouplingHomogen(dists, f, rhos, coilsType)
    for(j=KS; j<K; j++){
        const lamda = Math.exp(-Y[j])
        const lamda2 = lamda ** 2
        

        let [R, P] = Ro(rhos, f, lamda, thicks)
        
        
        if (coilsType == 'H'){
            const input = math.multiply(-1, math.pow(dists, 2), lamda2, math.subtract(R, P))
            
            Z = math.add(Z, math.multiply(input, ck.H[j]))
        }
        else{
            const input = math.multiply(-1, lamda, dists, R)
            Z = math.add(Z, math.multiply(input, ck.V[j]))
        }

        
        

    }
    let appConduc = 4*Z.im / (MU * 2 * Math.PI * f * dists ** 2)
    return(appConduc * 1000)
}


function div(w, z){
    const low = 1 / (z.re ** 2 + z.im ** 2)

    const re = w.re*z.re + w.im*z.im
    const im = w.im*z.re - w.re*z.im
    const upp = math.complex(re, im)

    const wz = math.multiply(low, upp)
    return wz
}


function row(){
    const cellTable = document.createElement('td')
    const inputTable = document.createElement('input')

    // Add CSS
    inputTable.setAttribute("style", "height: 20px; width:70px ; margin :0px")
    inputTable.type = 'number'
    inputTable.min = 0
    inputTable.max = 100000

    // Append
    cellTable.appendChild(inputTable)

    return cellTable

}

function addRowTable(numLayer){
    const bInputModel = document.getElementById('bInputModel')

    // creat elements
    const rowTable = document.createElement('tr')
    const headTable = document.createElement('th')
    const layer = document.createTextNode(numLayer)

    const cellTable1 = row()
    const cellTable2 = row()

    headTable.appendChild(layer)
    rowTable.appendChild(headTable)
    rowTable.appendChild(cellTable1)
    rowTable.appendChild(cellTable2)    
 
    bInputModel.appendChild(rowTable)

}


function removeRowTable(){
    const rowTable = document.querySelectorAll('#inputObs tr')
    const numLayer = rowTable.length - 1
    if (numLayer > 1){
        rowTable[numLayer].remove()
    }
    else{
        alert('One layer required! You cannot remove row!')
    } 
       
}

function getModel(){
    const rows = document.querySelectorAll('#inputObs input')
    let thickness = []
    let rhos = []
    rows.forEach((item, i) => {
        if (i % 2 === 0){
            thickness.push(parseFloat(item.value))
        }
        else{
            rhos.push(parseFloat(item.value))
        }
    });
    return [thickness, rhos]
}


function getDataCalc(){
    const rowsH = document.querySelectorAll('#inputH input')
    const rowsV = document.querySelectorAll('#inputV input')

    let inputDataCalH = []
    let inputDataCalV = []

    rowsH.forEach(e =>{
        inputDataCalH.push(parseFloat(e.value))
    })

    rowsV.forEach(e =>{
        inputDataCalV.push(parseFloat(e.value))
    })

    return [inputDataCalH, inputDataCalV]

}

function setCalculateData(data){
    const calcData = document.querySelectorAll("#calcData td")
    calcData.forEach((element, index) => {
        element.innerHTML = data[index].toPrecision(6)
    });
}


function validateInput(data){
    data.forEach(element =>{
        if (isNaN(element)) return 0
        if (element < 0) return 0
    })
}

const btnPlus = document.getElementById('btnPlus')
btnPlus.addEventListener('click', function(){
    numLayer = document.querySelectorAll('#inputObs tr').length
    if (numLayer <= 20){
        addRowTable(numLayer)
    }
    else{
        alert('Maximum adding row only 20!')
    }
    
})

const btnMinus = document.getElementById('btnMinus')
btnMinus.addEventListener('click', function(){
    removeRowTable()
})

const btnExecute = document.getElementById('btnExecute')
btnExecute.addEventListener('click', function(){

    const [thickness, rhos] = getModel()
    const [inputDataCalH, inputDataCalV] =  getDataCalc()

    const inputModel = thickness.concat(rhos)
    const inputDataObs = inputDataCalH.concat(inputDataCalV)

    let val = true
    inputModel.forEach(element =>{
        if(isNaN(element) || element < 0){
            val = false
        }
    })

    inputDataObs.forEach(element =>{
        if(isNaN(element) || element < 0){
            val = false
        }
    })

    if (val){
        const dists     = [10, 20, 40]
        const f         = [6400, 1600, 400]

        const H = new Array(dists.length).fill(0)
        const V = new Array(dists.length).fill(0)

        for (i=0; i<dists.length;i++){
            H[i] = condApp(f[i], thickness, rhos, dists[i], 'H')
            V[i] = condApp(f[i], thickness, rhos, dists[i], 'V')
        }

        const data = H.concat(V)
        setCalculateData(data)
        }
    else{
        alert('Input Error!! Empty or minus data input is not allowed!!')
    }
    
})


// set the dimensions and margins of the graph
var margin = {top: 10, right: 30, bottom: 30, left: 60},
    width = 460 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

//Read the data
d3.csv("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/2_TwoNum.csv", function(data) {

  // Add X axis
  var x = d3.scaleLinear()
    .domain([0, 4000])
    .range([ 0, width ]);
  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(x));

  // Add Y axis
  var y = d3.scaleLinear()
    .domain([0, 500000])
    .range([ height, 0]);
  svg.append("g")
    .call(d3.axisLeft(y));

  // Add dots
  svg.append('g')
    .selectAll("dot")
    .data(data)
    .enter()
    .append("circle")
      .attr("cx", function (d) { return x(d.GrLivArea); } )
      .attr("cy", function (d) { return y(d.SalePrice); } )
      .attr("r", 1.5)
      .style("fill", "#69b3a2")

})