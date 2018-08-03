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
        noise_type="determ"
        originalObsErr = obsErr = 5,
        freq = 200;
        std=2.0

   //global L63 parameters
    var rho=28.0;
        sigma=10.0;
        beta=8.0/3.0;


    // Truth state vector
    var truth_i = [0.0,0.0,0.0], truth_f, truth_i_plot;
    truth_i_plot = truth_i;

    window.addEventListener('resize', function(){two.clear()},false);
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
            UpdateFromUI()
            // Step truth forward
            if (noise_type=="determ"){
            truth_f = step(truth_i);
            }
            if (noise_type=="stoch"){
              noise=WhiteNoise()
              //adds white noise onto deterministic update
              truth_f =step(truth_i).map(function(num,idx){return num+std*Math.sqrt(dt)*noise[idx];});
            }
            if (noise_type=="hi_stoch"){
              //make noise in EOF space, zero eof 3, and map back to x-space
              noise=WhiteNoise()
              noise[2]=0
              noise=inverse_eof_transform(noise)
              //adds white noise onto deterministic update
              truth_f =step(truth_i).map(function(num,idx){return num+std*Math.sqrt(3*dt/2)*noise[idx];});
            }
            if (noise_type=="lo_stoch"){
              //make noise in EOF space, zero eofs 1-2, and map back to x-space
              noise=WhiteNoise()
              noise[0]=noise[1]=0
              noise=inverse_eof_transform(noise)
              //adds white noise onto deterministic update
              truth_f =step(truth_i).map(function(num,idx){return num+std*Math.sqrt(3*dt)*noise[idx];});
            }
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
    function othername() {
        var input = document.getElementById("userInput").value;
        alert(input)
    }
    function UpdateFromUI(){
        noise_type=doc.getElementById('noise_mode').value
      }

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
    // Runge-Kutta integration
    function step(prev) {
        var k1 = ode(prev);
        var k2 = ode(math.add(prev, math.multiply(0.5*dt,k1)));

        return math.add(prev, math.multiply(dt,k2));
    }

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
