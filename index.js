//html elements
//canvases (julia canvas (default) and the point canvas)
const canvas = document.getElementById('canvas')
const ctx = canvas.getContext('2d');

canvas.width = 720*1.5;
canvas.height = 360*1.5;

ctx.translate(0, canvas.height);
ctx.scale(1, -1);

const pcanvas = document.getElementById('pcanvas');
const pctx = pcanvas.getContext('2d');

pcanvas.width = 720*1.5;
pcanvas.height = 360*1.5;

pctx.translate(0, canvas.height);
pctx.scale(1, -1);


//variables
var points; //points to iterate
var nextpoints; //this is so calculations can be made ahead of time for smooth animation

//divs
var parameterDivs = document.getElementsByClassName('parameter');
var xdensediv = document.getElementById("xdensediv");
var ydensediv = document.getElementById("ydensediv");

//inputs
var xdenseinput = document.getElementById("xdense");
var ydenseinput = document.getElementById("ydense");


//sliders (remember to use slider.oninput = ()->{} to define what it does)



//helper functions

function to_zoom(value){
    let scale = Math.round(value)/100;
    return scale;
}

function rgbToHex(rgb){
    let x = Math.floor(rgb[0]/16);
    let x2 = parseInt(rgb[0]-x*16);
    let y = Math.floor(rgb[1]/16);
    let y2 = parseInt(rgb[1]-y*16);
    let z = Math.floor(rgb[2]/16);
    let z2 = parseInt(rgb[2]-z*16);
    x = x.toString(16)
    x2 = x2.toString(16)
    y = y.toString(16)
    y2 = y2.toString(16)
    z = z.toString(16)
    z2 = z2.toString(16)

    return '#' + x + x2 + y + y2 + z + z2;
}

function plotArray(canvas, array, strength=120){

    let ctx = canvas.getContext('2d');

    let xres = canvas.width / array.length;
    let yres = canvas.height / array[0].length;

    for(let i=0; i<canvas.width; i+=xres){
        for(let j=0; j<canvas.height; j+=yres){
            let color = array[parseInt(i/xres)][parseInt(j/yres)];
            let style = rgbToHex([color*strength, 0, 0]);
            ctx.fillStyle = style;
            ctx.fillRect(i, j, Math.ceil(xres)+1, Math.ceil(yres)+1);
            ctx.fill();
        }
    }
}

function plotPoint(canvas, z, z0, z1, s, color = 'gold'){

    let ctx = canvas.getContext('2d');

    let diff = z1.sub(z0);

    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc((z.re-z0.re)/diff.re * canvas.width, (z.im-z0.im)/diff.im * canvas.height, s, 0, 2*Math.PI);
    ctx.fill();
}

function scatter(canvas, array, z0, z1, s, color="gold"){
    let ctx = canvas.getContext('2d');
    ctx.fillStyle = color;

    for(let z of array){
        let diff = z1.sub(z0);

        ctx.beginPath();
        ctx.arc((z.re-z0.re)/diff.re * canvas.width, (z.im-z0.im)/diff.im * canvas.height, s, 0, 2*Math.PI);
        ctx.fill();
    }
}


//HTML PAGE FUNCTIONS

function biterateFunc(){

    if(interval){
        return;
    }

    let t = 0;
    let tmax = 2;

    let temp = points;

    interval = setInterval(()=>{
        
        pctx.clearRect(0, 0, pcanvas.width, pcanvas.height);

        if(t>=tmax){

            points = nextpoints;
            nextpoints = iterate(points);
            scatter(pcanvas, points, minz, maxz, s);
            interval = clearInterval(interval);

        }

        else{

            let ct = new Complex((t/tmax)**3, 0);
            let ct1 = new Complex(1-(t/tmax)**3, 0);

            
            scatter(pcanvas, temp, minz, maxz, s);

            for(let i=0; i<points.length; i++){
                temp[i] = points[i].mult(ct1).add(nextpoints[i].mult(ct));
            }


            t+=16/1000;
        }

    }, 16)

}

function juliaFunc(){

    // with polynomial original julia function
    // let R = (z)=>{return evaluate(P, z)};
    // J = polynomial_julia(R, eps2, minz, maxz, args = {'bound':2,'niter':40, 'halo':0})
    // var e = binaryEdge(J);
    // plotArray(canvas, e.data, strength = 120);


    //with rational julia function
    J = julia(P, Q, juliaN);

    for(let z of J){
        plotPoint(canvas, z, minz, maxz, 1, "red");
    }

    return
}

function pointsFunc(){
    let X = jsn.linspace(minz.re, maxz.re, xdense);
    let Y = jsn.linspace(minz.im, maxz.im+0.1, ydense); //+0.1 is arbitrary, just for prettiness
    points = [];

    for(let x of X){
        for(let y of Y){
            points.push(new Complex(x, y));
        }
    }

    scatter(pcanvas, points, minz, maxz, s);

    nextpoints = iterate(points);

    return
}


//update functions

function updateDensity(){

    xdense = parseInt(xdenseinput.value);
    ydense = parseInt(ydenseinput.value);

    xdensediv.innerHTML = "x points density: " + xdenseinput.value;
    ydensediv.innerHTML = "y points density: " + ydenseinput.value;

}

function updatePoly(){

}

function update(){
    
}
