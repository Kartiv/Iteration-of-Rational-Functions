//html elements

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

const biterate = document.getElementById('biterate');

const colorDict = {0:'black', 1:'lightblue', 2:'magenta', 3:'gold', 4:'#808000'}

//variables
let interval;

//functions

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

function julia(R, eps, z0, z1, args = {'bound':4, 'niter':20, 'halo':1}){

    let bound;
    let niter;
    let halo;

    if(args['bound']==undefined){
        bound = 4;
    }
    else{
        bound = args['bound'];
    }
    if(args['niter']==undefined){
        niter = 20;
    }
    else{
        niter = args['niter'];
    }
    if(args['halo']==undefined){
        halo = true;
    }
    else{
        halo = args['halo'];
    }

    let X = jsn.arange(z0.re, z1.re, eps);
    let Y = jsn.arange(z0.im, z1.im, eps);

    let C = [];

    for(let x of X){
        C.push([]);
        for(let y of Y){
            let flag = 0;
            let z0 = new Complex(x, y);

            for(let i=0; i<niter; i++){
                z0 = R(z0);
                if(z0.abs>bound*bound){
                    flag = halo*i+1;
                    break
                }
            }
            C[C.length-1].push(flag);
        }
    }
    return new Matrix(C);
}

function binaryEdge(picture){
    let C = picture.data;
    let edgeCoords = [];

    for(let i=1; i<C.length-1; i++){
        for(let j=1; j<C[0].length-1; j++){

            let flag = true;
            let sub = picture.subMatrix(i-1, i+1+1, j-1, j+1+1).data;
            for(let r=0; r<sub.length; r++){
                for(let c=0; c<sub[0].length; c++){
                    if(sub[r][c]!=C[i][j]){
                        flag = false;
                    }
                }
            }

            if(!flag){
                edgeCoords.push([i,j])
            }
        }
    }

    let im = [];
    for(let i=0; i<C.length; i++){
        im.push([]);
        for(let j=0; j<C[0].length; j++){
            im[i].push(0);
        }
    }
    for(let p of edgeCoords){
        im[p[0]][p[1]] = 1;
    }

    return new Matrix(im);
}

function iterate(points){
    let newpoints = [];
    for(let p of points){
        newpoints.push(R(p));
    }
    return newpoints;
}

function biterateFunc(){

    if(interval){
        return;
    }

    let newpoints = iterate(points);

    let temp = [];

    let t = 0;
    let tfinal = 1;

    interval = setInterval(()=>{

        pctx.clearRect(0,0,pcanvas.width, pcanvas.height);

        if(t>=1){
            for(let i=0; i<newpoints.length; i++){
                plotPoint(pcanvas, newpoints[i], z0, z1, s);
            }
            points = newpoints;
            interval = clearInterval(interval);
        }

        for(let i=0; i<newpoints.length; i++){
            temp[i] = points[i].mult(new Complex(1-t**3, 0)).add(newpoints[i].mult(new Complex(t**3, 0)));
            plotPoint(pcanvas, temp[i], z0, z1, s);
        }

        t+=16/(1000*tfinal);

    }, 16)
    for(let p of points){
        plotPoint(pcanvas, p, z0, z1, s);
    }
}

function R(z){
    return (z.mult(z).add((new Complex(0.33, 0.175)))); //remember 0.3,0.5 since it can be a good component tester
    //nice functions: c=0.3+0.5i, c=-0.1-0.7i, c=-0.5-0.51i
}


//plotting the julia set

let z0 = new Complex(-2, -1.2);
let z1 = new Complex(2, 1.2);
let eps = 0.002; //0.002 good

let C = julia(R, eps, z0, z1, args = {'bound':4,'niter':100, 'halo':0});

let e = binaryEdge(C);

plotArray(canvas, e.data, strength = 120);


//Showing the dynamics
 
let s = 1; //width of points

let xdense = 100; //how many points per row
let ydense = 50; //how many points per column

let X = jsn.linspace(z0.re, z1.re, xdense);
let Y = jsn.linspace(z0.im, z1.im+0.1, ydense); //+0.1 is arbitrary, just for prettiness
let points = [];

for(let x of X){
    for(let y of Y){
        points.push(new Complex(x, y));
    }
}

for(let p of points){
    plotPoint(pcanvas, p, z0, z1, s);
}

