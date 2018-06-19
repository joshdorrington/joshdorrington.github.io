(function(doc,win) {
    // Create two.js element
    var el = document.getElementById("main"),
        two = new Two({
            fullscreen: true
        }).appendTo(el);

    //two.renderer.domElement.style.backgroundColor = 'rgb(80,80,80)';

    // Global parameters
    var dt = 0.0001,
        nEns = 80,
        assim_freq = 55,
        originalObsErr = obsErr = 5,
        freq = 2000;

   //global CdV parameters
    var C=0.1,
        b=0.5,
        g=0.2,
        x1f=0.95,
        x4f=-0.76095,
        epsilon=16*Math.sqrt(2)/(5*Math.PI),
        beta=1.25;

    //and the derived spectral coefficients
    var beta1=beta*b*b/(1+b*b),
        beta2=beta*b*b/(4+b*b),
        alpha1=1*8*b*b*Math.sqrt(2)/(3*Math.PI*(1+b*b)),
        alpha2=4*8*(3+b*b)*Math.sqrt(2)/(15*Math.PI*(4+b*b)),
        gammatilde1=g*1*4*b*Math.sqrt(2)/(3*Math.PI),
        gammatilde2=g*2*4*b*Math.sqrt(2)/(15*Math.PI),
        gamma1=g*4*b*Math.sqrt(2)/(3*Math.PI*(1+b*b)),
        gamma2=g*32*b*Math.sqrt(2)/(15*Math.PI*(4+b*b)),
        delta1=64*Math.sqrt(2)*b*b/(15*Math.PI*(1+b*b)),
        delta2=64*Math.sqrt(2)*(b*b-3)/(15*Math.PI*(4+b*b));

    // Truth state vector
    var truth_i = [0.927796,0.160574,-0.0969953,-0.672758,-0.118184,0.214505], truth_f, truth_i_plot;
    truth_i_plot = truth_i;

    el.addEventListener('click', function(e) {
        console.log(e.clientX);
        console.log(e.clientY);
        truth_i = truth_i_plot = mapToCoord([e.clientX, e.clientY]);
        two.clear();
    }, false);

    var handles = [];

    // Main loop
    two.bind("update", function(frameCount) {
        for (var i = 0; i < freq; i++) {
            // Step truth forward
            truth_f = step(truth_i);
            truth_i = truth_f;
        }

        // Update truth
        handles.push(plotTruth(truth_i_plot, truth_f));

        // Update opacity
        handles = updateOpacity(handles);
        truth_i_plot = truth_i;
    });

    // Run main loop
    setInterval(function() {
        two.update();
    }, 24);

    function plotTruth(loc_i, loc_f) {
        loc_i_px = mapToPx(loc_i);
        loc_f_px = mapToPx(loc_f);
        var handle = two.makeLine(loc_i_px.x, loc_i_px.y, loc_f_px.x, loc_f_px.y);
        handle.stroke = 'rgb(255,255,255)';
        handle.linewidth = 2.0;
        handle.cap = 'round';
        return handle;
    }

    function updateOpacity(handles) {
        var numKeep = 1000;
        for (var i = 0; i < handles.length; i++) {
            handles[i].opacity -= 1/numKeep;
        }
        if (handles.length > numKeep) handles.shift();
        return handles;
    }

    function mapToPx(state) {

        //applies eof transform to data
        var eofs=math.matrix([[0.24480295,-0.16261957,0.54734502,-0.47431644,-0.0384,-0.6225666 ],
        [0.0496688 ,-0.44083483,-0.21379727,-0.25879109,0.82575797,0.09294802]])
        var projected_state=math.multiply(eofs,state)

        return {
            x: (projected_state._data[0]+0.28)*two.width/1.024,
            y: (0.685-projected_state._data[1])*two.height/1.023
        };
    }

    function mapToCoord(loc_px) {

        //applies inverse eof transform
        var eof_inverse=math.matrix([[ 0.24480295,0.0496688,-0.43292675,0.62904923,0.09100189,0.58838311],
        [-0.16261957, -0.44083483, -0.23104762, -0.49480442, -0.44425248,0.53258318],
        [ 0.54734502, -0.21379727, -0.26855881,  0.14640603, -0.57997197,-0.4741077 ],
        [-0.47431644, -0.25879109,  0.45143826,  0.56826685, -0.42571605,0.00965485],
        [-0.0384    ,  0.82575797,  0.05266652, -0.09479475, -0.52607519,0.16773285],
        [-0.6225666 ,  0.09294802, -0.69317914,  0.07821558, -0.00128128,-0.3422789 ]])
        //the 4 smallest eof modes are set to the climatological mean value
        var scaled_px=[(1.024*loc_px[0]/two.width)-0.28,(-1.023*loc_px[1]/two.height)+0.685,-0.7325,0.1736,0.3842,0.5871]
        return math.multiply(eof_inverse,scaled_px)._data
    }

    // Runge-Kutta integration
    function step(prev) {
        var k1 = ode(prev);
        var k2 = ode(math.add(prev, math.multiply(0.5*dt,k1)));

        return math.add(prev, math.multiply(dt,k2));
    }

    function ode(state) {
        var x1 = state[0],
            x2 = state[1],
            x3 = state[2],
            x4 = state[3],
            x5 = state[4],
            x6 = state[5];


        return [dX1dT(x1,x2,x3,x4,x5,x6), dX2dT(x1,x2,x3,x4,x5,x6), dX3dT(x1,x2,x3,x4,x5,x6),dX4dT(x1,x2,x3,x4,x5,x6),dX5dT(x1,x2,x3,x4,x5,x6),dX6dT(x1,x2,x3,x4,x5,x6)];
    }

    // Lorenz '63 equations
    function dX1dT(x1,x2,x3,x4,x5,x6) {
        return -C*(x1-x1f)+gammatilde1*x3
    }

    function dX2dT(x1,x2,x3,x4,x5,x6){

        return -C*(x2)+beta1*x3-alpha1*x1*x3-delta1*x4*x6
    }

    function dX3dT(x1,x2,x3,x4,x5,x6){

        return -C*(x3)-beta1*x2-gamma1*x1+alpha1*x1*x2+delta1*x4*x5
    }
     function dX4dT(x1,x2,x3,x4,x5,x6){
        var σ = 10.0;
        return -C*(x4-x4f)+gammatilde2*x6+epsilon*(x2*x6-x3*x5)
    }
     function dX5dT(x1,x2,x3,x4,x5,x6){
        var σ = 10.0;
        return -C*(x5)+beta2*x6-alpha2*x1*x6-delta2*x3*x4
    }
     function dX6dT(x1,x2,x3,x4,x5,x6){
        var σ = 10.0;
        return -C*(x6)-beta2*x5-gamma2*x4+alpha2*x1*x5+delta2*x2*x4
    }

}(document, window));
