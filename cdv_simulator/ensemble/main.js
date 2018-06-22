(function(doc,win) {
    // Create two.js element
    var el = document.getElementById("main"),
        two = new Two({
            fullscreen: true
        }).appendTo(el);

    //two.renderer.domElement.style.backgroundColor = 'rgb(29,29,29)';

    // Global parameters
    var dt = 0.005,
        nEns = 200,
        noise_type="determ"
        originalObsErr = obsErr = 5,
        freq = 20;

   //global CdV parameters
    var C=0.1,
        b=0.5,
        g=0.2,
        x1f=0.95,
        x4f=-0.76095,
        epsilon=16*Math.sqrt(2)/(5*Math.PI),
        beta=1.25,
        sigma=0.007;

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

        // Flags
        var clicked = false;

    // Initialise array storing ensemble of state vectors
    var ensemble_i = [], ensemble_f;
    for (var i = 0; i < nEns; i++) {
        ensemble_i.push([0,0,0,0,0,0]);
    }

    // Initialise SVG objects
    var ensembleHandle = [];
    for (var i = 0; i < nEns; i++) {
        ensembleHandle.push(two.makeCircle(0.0,0.0,8.0, 1.0));
        ensembleHandle[i].fill = 'rgba(0,0,0,0.3)';
        ensembleHandle[i].linewidth = 0;
    }

    // Event listener for clicks
    el.addEventListener('click', function(e) {
        if (!clicked) {

            // Main loop
            two.bind("update", function(frameCount) {
                // Step ensemble forward
                for (var t=0; t<freq; t++){
                  UpdateFromUI()
                  for (var i = 0; i < nEns; i++) {

                    if (noise_type=="determ"){
                    ensemble_f[i] = step(ensemble_i[i]);
                    }
                    if (noise_type=="stoch"){
                      noise=WhiteNoise()
                      //adds white noise onto deterministic update
                      ensemble_f[i] = step(ensemble_i[i]).map(function(num,idx){return num+sigma*Math.sqrt(dt)*noise[idx];});
                    }
                    if (noise_type=="hi_stoch"){
                      //make noise in EOF space, zero eofs 4-6, and map back to x-space
                      noise=WhiteNoise()
                      noise[3]=noise[4]=noise[5]=0
                      noise=inverse_eof_transform(noise)
                      //adds white noise onto deterministic update
                      ensemble_f[i] = step(ensemble_i[i]).map(function(num,idx){return num+sigma*Math.sqrt(dt)*noise[idx];});
                    }
                    if (noise_type=="lo_stoch"){
                      //make noise in EOF space, zero eofs 1-3, and map back to x-space
                      noise=WhiteNoise()
                      noise[0]=noise[1]=noise[2]=0
                      noise=inverse_eof_transform(noise)
                      //adds white noise onto deterministic update
                      ensemble_f[i] = step(ensemble_i[i]).map(function(num,idx){return num+sigma*Math.sqrt(dt)*noise[idx];});
                    }
                  }
                if (t==freq-1){updateEnsemble(ensembleHandle, ensemble_i, ensemble_f);}
                ensemble_i = ensemble_f.slice()
                }



            });

            // Make ensemble visible
            for (var i = 0; i < nEns; i++) {
                ensembleHandle[i].fill = 'rgba(255,255,255,0.3)';
            }
            clicked = true;
        }

        var centre = mapToCoord([e.clientX, e.clientY]);

        // Initialise ensemble members as a random spherical cloud surrounding click location
        for (var i = 0; i < nEns; i++) {
            var u = Math.random(), v = Math.random();
            var theta = 2*Math.PI*u;
            ensemble_i[i][0] = centre[0] + (Math.random()-0.5)/100
            ensemble_i[i][1] = centre[1] + (Math.random()-0.5)/100
            ensemble_i[i][2] = centre[2] + (Math.random()-0.5)/100
            ensemble_i[i][3] = centre[3] + (Math.random()-0.5)/100
            ensemble_i[i][4] = centre[4] + (Math.random()-0.5)/100
            ensemble_i[i][5] = centre[5] + (Math.random()-0.5)/100
        }
        ensemble_f = ensemble_i.slice();
    }, false);

    // Run main loop
    setInterval(function() {
       two.update();
    }, 24);

    //box muller, returns normally distributed random number
    function randn_bm() {
        var u = 0, v = 0;
        while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
        while(v === 0) v = Math.random();
        var R=Math.sqrt( -2.0 * Math.log( u ))
        return  [R*Math.cos(2.0*Math.PI*v),R*Math.sin(2.0*Math.PI*v)];
    }
    function WhiteNoise(){
       var noise_vec = new Array(6)
       for (var i = 0; i < 2+(noise_vec.length)/2; i=i+2) {
           arr=randn_bm()
           noise_vec[i]=arr[0]
           noise_vec[i+1]=arr[1]
       }
       return noise_vec
    }
    function UpdateFromUI(){
        noise_type=doc.getElementById('noise_mode').value
      }
    // Update locations of ensemble members
    function updateEnsemble(handles, locs_i, locs_f) {
        for (var i = 0; i < nEns; i++) {
            var loc_px = mapToPx(locs_f[i])
            handles[i].translation.set(loc_px.x, loc_px.y);

        }
    }

    // Map a state vector to pixels
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

    // Map position in pixels to state vector
    function mapToCoord(loc_px) {
        //the 4 smallest eof modes are set to the climatological mean value
        var scaled_px=[(1.024*loc_px[0]/two.width)-0.28,(-1.023*loc_px[1]/two.height)+0.685,-0.7325,0.1736,0.3842,0.5871]
        return inverse_eof_transform(scaled_px)
    }

    //applies inverse eof transform
    function inverse_eof_transform(input_vector){
      var eof_inverse=math.matrix([[ 0.24480295,0.0496688,-0.43292675,0.62904923,0.09100189,0.58838311],
      [-0.16261957, -0.44083483, -0.23104762, -0.49480442, -0.44425248,0.53258318],
      [ 0.54734502, -0.21379727, -0.26855881,  0.14640603, -0.57997197,-0.4741077 ],
      [-0.47431644, -0.25879109,  0.45143826,  0.56826685, -0.42571605,0.00965485],
      [-0.0384    ,  0.82575797,  0.05266652, -0.09479475, -0.52607519,0.16773285],
      [-0.6225666 ,  0.09294802, -0.69317914,  0.07821558, -0.00128128,-0.3422789 ]])

      return math.multiply(eof_inverse,input_vector)._data
    }

    // Compute the orientation of the difference vector between the two given vectors
    function getRotation(locs_i, locs_f) {
        var diffX = locs_f[0] - locs_i[0];
        var diffY = locs_f[1] - locs_i[1];

        return -Math.atan(diffY/diffX);
    }

    // Runge-Kutta integration
    function step(prev) {
        var k1 = ode(prev);
        var k2 = ode(math.add(prev, math.multiply(0.5*dt,k1)));

        return math.add(prev, math.multiply(dt,k2));
    }

    // The ODE for CdV model
    function ode(state) {
        var x1 = state[0],
            x2 = state[1],
            x3 = state[2],
            x4 = state[3],
            x5 = state[4],
            x6 = state[5];


        return [dX1dT(x1,x2,x3,x4,x5,x6), dX2dT(x1,x2,x3,x4,x5,x6), dX3dT(x1,x2,x3,x4,x5,x6),dX4dT(x1,x2,x3,x4,x5,x6),dX5dT(x1,x2,x3,x4,x5,x6),dX6dT(x1,x2,x3,x4,x5,x6)];
    }

    // CdV equations
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
        return -C*(x4-x4f)+gammatilde2*x6+epsilon*(x2*x6-x3*x5)
    }
     function dX5dT(x1,x2,x3,x4,x5,x6){
        return -C*(x5)+beta2*x6-alpha2*x1*x6-delta2*x3*x4
    }
     function dX6dT(x1,x2,x3,x4,x5,x6){
        return -C*(x6)-beta2*x5-gamma2*x4+alpha2*x1*x5+delta2*x2*x4
    }
}(document, window));
