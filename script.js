//variables
let interval;

var J;

let minz = new Complex(-2, -2);
let maxz = new Complex(2, 2);
let eps = 0.05; //epsilon for normal julia, honestly i have no idea why this epsilon matters so much for computation
let eps2 = 0.002; //epsilon for polynomial_julia

var xdense = 100;
var ydense = 50;

var s = 1; //size of points

var juliaN = 70000; //how many preimages to calculate
var newtonN1 = 40; //how many times to iterate newton raphson when calculating preimages
var newtonN2 = 100; //how many times to iterate newton raphson when calculating fixed points

//functions

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

function polynomial_julia(R, eps, z0, z1, args = {'bound':4, 'niter':20, 'halo':1}){

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

//Newton Raphson Helper Functions

function createPoly(arr){//turns nx2 array into degree n complex polynomial

    let P = [];
    for(let coeff of arr){
        P.push(new Complex(coeff[0], coeff[1]));
    }

    return P;
}

function evaluate(P, z){ //evaluates polynomial P at z

    let s = P[0];
    for(let i=1; i<P.length; i++){
        s=s.add(P[i].mult(z.pow(i)));
    }

    return s;
}

function derivative(P){ //returns derivative of polynomial

    if(P.length==1){
        return [new Complex(0,0)];
    }

    let dP = [];
    for(let i=1; i<P.length; i++){
        dP.push(P[i].mult(new Complex(i, 0)));
    }

    return dP;

}

function newton_raphson(P, z0, niter){
    let z = z0;
    let dP = derivative(P);
    let q = rational(P, dP);

    for(let i=0; i<niter; i++){
        z = z.sub(q(z));
    }

    return z;
}

function rational(P, Q){ //creates rational function from polynomials, represented as coefficient arrays
    
    return (z)=>{

        if(z.re == Infinity){
            if(P.length>Q.length){ //this line here assumes that the polynomials are represented without trailing zeros
                return z;
            }
            else if(P.length<Q.length){
                return new Complex(0,0);
            }
            else{
                return P[P.length-1].div(Q[Q.length-1]);
            }
        }

        else{

            let pz = evaluate(P, z);
            let qz = evaluate(Q, z);

            if(qz.abs<0.000000001){
                return new Complex(Infinity, Infinity);
            }

            else{
                return pz.div(qz);
            }
        }
    }

}

function remove_trailing_zeros(P){

    let n = 0;
    for(let i=P.length-1; i>0; i--){
        if(P[i].abs > 0.000000001){
            n = i;
            break;
        }
    }

    return jsn.parse(P, 0, n+1, 1);
}

function same_length(P,Q){
    
    let n = P.length;
    let m = Q.length;

    if(n<m){
        for(let i=0; i<m-n; i++){
            P.push(new Complex(0,0));
        }
    }
    else if(n>m){
        for(let i=0; i<n-m; i++){
            Q.push(new Complex(0,0));
        }
    }
}

function multiply_poly(P,Q){ //multiply P and Q

    same_length(P, Q);

    let vP = jsn.fft(jsn.fpad(P, 1));
    let vQ = jsn.fft(jsn.fpad(Q, 1));

    let prod = [];
    for(let i=0; i<vP.length; i++){
        prod.push(vP[i].mult(vQ[i]));
    }

    return remove_trailing_zeros(jsn.ifft(prod));
}

function subtract_poly(P, Q){

    same_length(P, Q);
    let h = [];

    for(let i=0; i<P.length; i++){
        h.push(P[i].sub(Q[i]));
    }

    return h;
}

function divide_poly(P, Q){ //assuming Q is a degree 1 polynomial which divides P

    same_length(P, Q);

    let vP = jsn.fft(jsn.fpad(P));
    let vQ = jsn.fft(jsn.fpad(Q));

    //locate zero of Q
    let n = -1;
    for(let i=0; i<vQ.length; i++){
        if(vQ[i].abs < 0.0000000000000001){
            n = i;
            break;
        }
    }

    //divide values

    let h = [];
    for(let i=0; i<vP.length; i++){
        if(i==n){
            h.push((vP[i].add(new Complex(0.00001, 0.00001))).div(vQ[i].add(new Complex(0.00001, 0.00001))));
        }
        else{
            h.push(vP[i].div(vQ[i]));
        }
    }

    return remove_trailing_zeros(jsn.ifft(h));

}

function julia_fixed_point(P, Q){ //finds fixed point of a rational function P/Q with multiplier 1 or with abs>1

    //calculate H = P-zQ
    let H = [P[0]];

    same_length(P, Q);

    for(let i=1; i<P.length; i++){
        H.push(P[i].sub(Q[i-1]));
    }

    H.push(Q[Q.length-1].mult(new Complex(-1, 0)));

    //calculate the derivative of R
    let dP = derivative(P);
    let dQ = derivative(Q);

    //make derivatives be same length as P, Q
    dP.push(new Complex(0,0));
    dQ.push(new Complex(0,0));

    let Qsquared = remove_trailing_zeros(multiply_poly(Q, Q));
    let numerator = remove_trailing_zeros(subtract_poly(multiply_poly(dP, Q), multiply_poly(dQ, P)));

    let dR = rational(numerator, Qsquared);
    
    //start looking for fixed points
    for(let i=0; i<P.length; i++){
        
        //generate guess for newton raphson method
        let x = jsn.random(-1, 1);
        let y = jsn.random(-1, 1);
    
        let z0 = new Complex(x,y);
        

        //calculate root of H
        let r1 = newton_raphson(H, z0, newtonN2);

        let v = dR(r1);
        
        if(v.abs>0 || Math.abs(v.abs-1)<0.00000000000001){
            return r1;
        }

        else{
            H = divide_poly(H, [r1, new Complex(1,0)]);
        }

    }

    return new Complex(Infinity, Infinity);
    
}

function julia(P, Q, niter){

    //Step 1: Calculate points in julia set, specifically itll be a fixed point that is repelling/with multiplier 1
    let z0 = julia_fixed_point(P, Q);

    let J = [z0];

    //step 2: Calculate a random pre-image of z0
    //note that by now P and Q are of the same length already
    for(let iter=0; iter<niter; iter++){

        //case where z0 is infinity - we find a root of Q
        if(z0.re == Infinity){
            //initialize newton-raphson guess
            let x = jsn.random(-5, 5);
            let y = jsn.random(-5 ,5);
        
            let w = new Complex(x,y);

            let z1 = newton_raphson(Q, w, newtonN1);
            J.push(z1);
            z0 = z1;
        }

        else{
            let h = [];
            for(let i=0; i<P.length; i++){
                h.push(P[i].sub((z0.mult(Q[i]))));
            }
        
            let x = jsn.random(-5, 5);
            let y = jsn.random(-5 ,5);
        
            let w = new Complex(x,y);
            let z1 = newton_raphson(h, w, newtonN1);
            J.push(z1);
            z0 = z1;
        }

    }

    return J;
}


// let P = createPoly([[1, 0], [0, 0], [1, 0]]);

//interesting functions:

//julia set is on a differentiable curve?
// let P = createPoly([[1, 0], [1, 0], [0, 1], [0,0]]);
// let Q = createPoly([[1, 0], [0, 1], [0, 1]]);

//the fatou cycles arent part of (exclusively) the biggest components
// let P = createPoly([[-1, 0], [1, 0], [2, 0], [0,0]]);
// let Q = createPoly([[1, 0], [0, 0], [0, 0]]);

//cycle i->infty->i , it just looks cool (also, indifferent fixed point it seems (siegel disc?))
// let P = createPoly([[-1, 0], [0, 0], [0, 1], [0,0]]);
// let Q = createPoly([[1, 0], [0, 0], [1, 0]]);

// //soccer ball
// let P = createPoly([[-1, 0], [0, 0], [1, 0], [0, 0]]);
// let Q = createPoly([[1, 0], [0, 0], [1, 0]]);

// //cool weird thing
// let P = createPoly([[-1, 0], [0, 0], [1, 0], [0, 0]]);
// let Q = createPoly([[1, 0], [0, 0], [0, 1]]);

let lamda = new Complex(Math.cos(2), Math.sin(2));
let alpha = new Complex(0.2, -0.21);

let P = [new Complex(0,0), new Complex(0,0), lamda, lamda.mult(alpha.conj())];
let Q = [alpha, new Complex(1,0)];



// let Q = createPoly([[1, 0]]);
let R = rational(P, Q);