(function(doc,win) {
    // Create two.js element
    var el = document.getElementById("main"),
        two = new Two({
            fullscreen: true
        }).appendTo(el);

    //two.renderer.domElement.style.backgroundColor = 'rgb(29,29,29)';

    // Global parameters
    var dt = 0.0001,
        nEns = 150,
        noise_type="determ"
        originalObsErr = obsErr = 5,
        freq = 200;
        std=3.0


    //global L63 parameters
     var rho=28.0;
         sigma=10.0;
         beta=8.0/3.0;

        // Flags
        var clicked = false;

    // Initialise array storing ensemble of state vectors
    var ensemble_i = [], ensemble_f,reg_swaps=[];
    for (var i = 0; i < nEns; i++) {
        ensemble_i.push([0,0,0]);
        reg_swaps.push(0)
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
                      ensemble_f[i] = step(ensemble_i[i]).map(function(num,idx){return num+std*Math.sqrt(dt)*noise[idx];});
                    }
                    if (noise_type=="hi_stoch"){
                      //make noise in EOF space, zero eof 3, and map back to x-space
                      noise=WhiteNoise()
                      noise[2]=0
                      noise=inverse_eof_transform(noise)
                      //adds white noise onto deterministic update
                      ensemble_f[i] = step(ensemble_i[i]).map(function(num,idx){return num+std*Math.sqrt(3*dt/2)*noise[idx];});
                    }
                    if (noise_type=="lo_stoch"){
                      //make noise in EOF space, zero eofs 1-2, and map back to x-space
                      noise=WhiteNoise()
                      noise[0]=noise[1]=0
                      noise=inverse_eof_transform(noise)
                      //adds white noise onto deterministic update
                      ensemble_f[i] = step(ensemble_i[i]).map(function(num,idx){return num+std*Math.sqrt(3*dt)*noise[idx];});
                    }
                    //change colour after regime switch
                    if ((reg_swaps[i]==0)&((ensemble_i[i][0]/ensemble_f[i][0])<0)){
                      reg_swaps[i]=1
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
            reg_swaps[i]=0
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
       var noise_vec = new Array(3)
       for (var i = 0; i < noise_vec.length; i=i+1) {
       		 arr=randn_bm()
           noise_vec[i]=arr[0]
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
            if (reg_swaps[i]==0){
              handles[i].fill='rgba(255,255,255,0.3)'
            }
            if (reg_swaps[i]==1){
              handles[i].fill='rgba(255,0,0,0.3)'
            }
        }
    }

    // Map a state vector to pixels
    function mapToPx(state) {

      return {
          x: ((state[0]+23)/46)*two.width,
          y: ((state[2]-56)/-57)*two.height
      };
    }

    //y is set to x on click
    function mapToCoord(loc_px) {
        var scaled_px=[(46*loc_px[0]/two.width)-23,(46*loc_px[0]/two.width)-23,(-57*loc_px[1]/two.height)+56]
        console.log(scaled_px)
        return scaled_px
    }

    //applies inverse eof transform
    function inverse_eof_transform(input_vector){
      var eof_inverse=math.matrix([[-0.66415,-0.00382,-0.74759],
      [-0.74759, -0.00895, 0.66416],
      [-0.00470,0.99999,-0.00093]])

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
      var x = state[0],
          y = state[1],
          z = state[2];

      return [dX1dT(x,y,z),dX2dT(x,y,z),dX3dT(x,y,z)];
  }

  // L63 equations
  function dX1dT(x,y,z) {
      return sigma*(y-x)
  }

  function dX2dT(x,y,z){

      return x*(rho-z)-y
  }

  function dX3dT(x,y,z){

      return x*y-beta*z
  }

}(document, window));
