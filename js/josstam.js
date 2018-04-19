
function josStamFluid(canvas) {
var iterations = 10;
var visc = .001;
var dt = 0.1;
var width;
var height;
var dim;
var grid_size;
var user_update;
var density;
var density_prev;
var velx;
var velx_prev;
var vely;
var vely_prev;
var button_id
var display;

this.get_id=function(id)
{
button_id = id;
}

//addition of force
function addForce(vel, vel0, dt)
{
	for (var i=0; i<grid_size ; i++ ) 
	{
		vel[i] += dt*vel0[i];
	}
}

//stable solver for density
function stable_solve(bound_num, current, previous, viscosity, divisor)
{
	if (viscosity === 0 && divisor === 1) 
	{
		for (var j=1 ; j<=height; j++)
		 {
			var row_num = j * dim;
			++row_num;
			for (var i = 0; i < width; i++)
			 {
				current[row_num] = previous[row_num];
				++row_num;
			 }
		}
		boundary(bound_num,current );
	} 
	else
	{
		var div = 1 / divisor;
		for (var num=0 ; num<iterations; num++) 
		{
			for (var j=1 ; j<=height; j++)
			 {
				var row_up = (j - 1) * dim;
				var row = j * dim;
				var row_down = (j + 1) * dim;
				var vel_current = current[row];
				++row;
				for (var i=1; i<=width; i++)
				{
					vel_current = current[row] = (previous[row] + viscosity*(vel_current+current[++row]+current[++row_up]+current[++row_down])) * div;
				}
			 }
			boundary(bound_num, current);
		}
	}
}

//diffuse density step
function diffuse_density(bound_num, current, previous, dt)
{
	var viscosity = visc*dt;
	var divisor=1 + 4*viscosity;
	stable_solve(bound_num, current, previous, viscosity, divisor);
}


//stable solver for velocity
function stable_solve_vel(curr_x, prev_x, curr_y, prev_y, viscosity, divisor)
{
	if (viscosity === 0 && divisor === 1) {
		for (var j=1 ; j <= height; j++) {
			var row_num = j * dim;
			++row_num;
			for (var i = 0; i < width; i++) {
				curr_x[row_num] = prev_x[row_num];
				curr_y[row_num] = prev_y[row_num];
				++row_num;
			}
		}
		boundary(1, curr_x);
		boundary(2, curr_y);
	} 
	else 
	{
		var div = 1/divisor;
		for (var itr=0 ; itr<iterations; itr++)
		 {
			for (var j=1 ; j <= height; j++) 
			{
				var row_up = (j - 1) * dim;
				var row = j * dim;
				var row_down = (j + 1) * dim;
				var vel_current_x = curr_x[row];
				var vel_current_y = curr_y[row];
				++row;
				for (var i = 1; i <= width; i++) {
					vel_current_x = curr_x[row] = (prev_x[row] + viscosity * (vel_current_x + curr_x[row] + curr_x[row_up] + curr_x[row_down])) * div;
					vel_current_y = curr_y[row] = (prev_y[row] + viscosity * (vel_current_y + curr_y[++row] + curr_y[++row_up] + curr_y[++row_down])) * div;
				}
			}
			boundary(1, curr_x);
			boundary(2, curr_y);
		}
	}
}


//diffusing velocity
function diffuse_vel(curr_x, prev_x, curr_y, prev_y, dt)
{
	var viscosity = visc*dt;
	stable_solve_vel(curr_x, prev_x, curr_y, prev_y, viscosity, 1 + 4 * viscosity);
}


//advect step
function advect(bound, current, previous, velx, vely, dt)
{
	for (var j = 1; j<= height; j++) 
	{
		var position = j * dim;
		for (var i = 1; i <= width; i++) 
		{
			++position;
			var positionx = i - dt * width * velx[position]; 
			var positiony = j - dt * height * vely[position];
			if (positionx < 1)
			{
				positionx = 1;
			}
			else if (positionx > width)
			{
				positionx = width;
			}
			  
			if (positiony < 1)
			{
				positiony = 1;
			}
			else if (positiony > height)
			{
				positiony = height;
			}
			
			var posx1= positionx | 0;
			var posx2 = posx1 + 1;
			var posy1 = positiony | 0;
			var posy2 = posy1 + 1;
			var interpx2 = positionx - Math.floor(posx1);
			var interpx1 = 1 - interpx2;
			var interpy2 = positiony - Math.floor(posy1);
			var interpy1 = 1 - interpy2;
			var row1 = posy1 * dim;
			var row2 = posy2 * dim;
			
			current[position] = interpx1 * (interpy1 * previous[posx1 + row1] + interpy2 * previous[posx1 + row2]) + interpx2 * (interpy1 * previous[posx2 + row1] + interpy2 * previous[posx2 + row2]);
		}
	}
	boundary(bound, current);
}

//project step
function project(curr_x, curr_y, prev_x, prev_y)
{
	var gradient = -0.5 / Math.sqrt(width * height);
	for (var j = 1 ; j <= height; j++ )
	{
		var row = j * dim;
		var row_up = (j - 1) * dim;
		var row_left = row - 1;
		var curr_row = row;
		var row_right = row + 1;
		var row_down = (j + 1) * dim;
		for (var i = 1; i <= width; i++ ) 
		{
			prev_y[++curr_row] = gradient * (curr_x[++row_right] - curr_x[++row_left] + curr_y[++row_down] - curr_y[++row_up]);
			prev_x[curr_row] = 0;
		}
	}
	boundary(0, prev_y);
	boundary(0, prev_x);
	stable_solve(0, prev_x, prev_y, 1, 4 );

	for (var j = 1; j<= height; j++ ) 
	{
		var curr = j * dim;
		var prev = j * dim - 1;
		var next = j * dim + 1;
		var row_up = (j - 1) * dim;
		var row_down = (j + 1) * dim;
		for (var i = 1; i<= width; i++) 
		{
			curr_x[++curr] -= 0.5 * width * (prev_x[++next] - prev_x[++prev]);
			curr_y[curr]   -= 0.5 * height * (prev_x[++row_down] - prev_x[++row_up]);
		}
	}
	boundary(1, curr_x);
	boundary(2, curr_y);
}

//density step 
function density_step(vel, vel0, vel_x, vel_y, dt)
{
	addForce(vel, vel0, dt);
	diffuse_density(0, vel0, vel, dt );
	advect(0, vel, vel0, vel_x, vel_y, dt );
}

//velocity step
function vel_step(curr_x, curr_y, prev_x, prev_y, dt)
{
	addForce(curr_x, prev_x, dt );
	addForce(curr_y, prev_y, dt );
	var temp = prev_x; 
	prev_x = curr_x; 
	curr_x = temp;
	var temp = prev_y; 
	prev_y = curr_y; 
	curr_y = temp;
	diffuse_vel(curr_x,prev_x,curr_y,prev_y, dt);
	project(curr_x, curr_y, prev_x, prev_y);
	var temp = prev_x; 
	prev_x = curr_x; 
	curr_x = temp; 
	var temp = prev_y; 
	prev_y = curr_y; 
	curr_y = temp;
	advect(1, curr_x, prev_x, prev_x, prev_y, dt);// with respect to x
	advect(2, curr_y, prev_y, prev_x, prev_y, dt);// with respect to y
	project(curr_x, curr_y, prev_x, prev_y );
}


//setting and getting fluid properties
function Fluid_prop(density, xvel, yvel) 
{

	this.setDensity = function(x, y, d) 
	{
		 density[(x ) + (y ) * dim] = d;
	}
	this.getDensity = function(x, y)
	{
		 return density[(x ) + (y ) * dim];
	}
	this.setVelocity = function(x, y, xv, yv) 
	{
		 xvel[(x ) + (y ) * dim] = xv;
		 yvel[(x ) + (y ) * dim] = yv;
	}
	this.getXVelocity = function(x, y) 
	{
		 return xvel[(x ) + (y ) * dim];
	}
	this.getYVelocity = function(x, y) 
	{
		 return yvel[(x) + (y ) * dim];
	}
	
	 this.clearCanvas=function()
	{
	dim = width + 2;
	grid_size = (width+2)*(height+2);
	density = new Array(grid_size);
	density_prev = new Array(grid_size);
	velx = new Array(grid_size);
	velx_prev = new Array(grid_size);
	vely = new Array(grid_size);
	vely_prev = new Array(grid_size);
	for (var i = 0; i < grid_size; i++)
	{
		density_prev[i] = velx_prev[i] = vely_prev[i] = density[i] = velx[i] = vely[i] = 0;
	  }
	  return;
	 }
	
}



//getting the values from user interaction
function getfrom_user(density_prev, velx_prev, vely_prev)
{
	for (var i = 0; i < grid_size; i++)
	{
		velx_prev[i] = vely_prev[i] = density_prev[i] = 0.0;
	}
	user_update(new Fluid_prop(density_prev, velx_prev, vely_prev));
}


// main simulation algorithm
this.update = function ()
{
	getfrom_user(density_prev, velx_prev, vely_prev);
	vel_step(velx, vely, velx_prev, vely_prev, dt);
	density_step(density, density_prev, velx, vely, dt);
	display(new Fluid_prop(density, velx, vely));
}


//updates the display
this.update_disp= function(new_disp) 
{
	display = new_disp; //velocity or density
}



this.set_user_update = function(result_user) 
{
	user_update = result_user;
}

//collision detection with boundary
function boundary(boundary_num, vel)
{
	if (boundary_num===1) 
	{
		for (var i = 1; i <= width; i++) 
		{
			var down_coord= i + (height+1) *dim;
			var down_coord1= i + height * dim;
			vel[i] =  vel[i + dim];
			vel[down_coord] = vel[down_coord1];
		}

		for (var j = 1; i <= height; i++)
		{
			var left=j * dim;
			var left_next=1 + j * dim;
			vel[left] = -vel[left_next];
			var right=(width + 1) + j * dim;
			var right_next=width + j * dim;
			vel[right] = -vel[right_next];
		}
	} 
	else if (boundary_num === 2)
	{
		for (var i = 1; i <= width; i++) 
		{
			var row_up=i + dim;
			var row_down= i + height * dim;
			var row_down1= i + (height + 1) * dim;
			vel[i] = -vel[row_up];
			vel[row_down1] = -vel[row_down];
		}

		for (var j = 1; j <= height; j++)
		{
			var col1=j * dim;
			var col2=1 + j * dim;
			var col3=(width + 1) + j * dim;
			var col4=width + j * dim;
			vel[col1] =  vel[col2];
			vel[col3] =  vel[col4];
		}
	} 
	else 
	{
		for (var i = 1; i <= width; i++) 
		{
			var row1 =i + dim;
			var row2 =i + (height + 1) * dim;
			var row3= i + height * dim;
			vel[i] =  vel[row1];
			vel[row2] = vel[row3];
		}

		for (var j = 1; j <= height; j++)
		{
			var col_new1=j * dim;
			var col_new2=1 + j * dim;
			var col_new3=(width + 1) + j * dim;
			var clo_new4=width + j * dim
			vel[col_new1]=vel[col_new2];
			vel[col_new3]=vel[clo_new4];
		}
	}
	var bottom_left = (height + 1) * dim;
	var topRight=(width+1);
	var bottomRight = (width+1)+bottom_left;
	var bottomRight_up = (width + 1) + height * dim;
	vel[0]= 0.5 * (vel[1] + vel[dim]);
	vel[bottom_left]= 0.5 * (vel[1 + bottom_left] + vel[height * dim]);
	vel[topRight] = 0.5 * (vel[width] + vel[topRight+dim]);
	vel[bottomRight] = 0.5 * (vel[width + bottom_left] + vel[bottomRight_up]);
}

//clearing the canvas
function clearCanvas()
{
	dim = width + 2;
	grid_size = (width+2)*(height+2);
	density = new Array(grid_size);
	density_prev = new Array(grid_size);
	velx = new Array(grid_size);
	velx_prev = new Array(grid_size);
	vely = new Array(grid_size);
	vely_prev = new Array(grid_size);
	for (var i = 0; i < grid_size; i++)
	{
		density_prev[i] = velx_prev[i] = vely_prev[i] = density[i] = velx[i] = vely[i] = 0;
	}
	frames=1;
	actual_frames=0;
}
this.clearCanvas = clearCanvas;
width = 128;
height = 128;
clearCanvas();
		
}


//rendering part starts here

(function () 
{
var secondCanvas;
var secondCanvasValues;
var canvas;

function makeSecondCanvas(fluid)
 {
	canvas =document.getElementById("canvas");
	canvas.style.background = 'white';
	secondCanvas = document.createElement("canvas");
	secondCanvas.width = 128;
	secondCanvas.height = 128;
	var context = secondCanvas.getContext("2d");
	try
	{
		secondCanvasValues = context.createImageData(secondCanvas.width, secondCanvas.height);
	} 
	catch(e) 
	{
		return null;
	}
   
	var data_length = secondCanvas.width * secondCanvas.height * 4;
	for (var i=3; i<data_length; i+=4)
	{
		secondCanvasValues.data[i] = 255;
   }
		 
	secondCanvasValues.data[0] = 0;
}

//cardinal color
function cardinal(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					else if(frames>=400 && frames<=1000)
					{
					data[ red_index] =  fluid.getDensity(x, y) * 255 / 6;
					if(frames==1000)
					{
					frames=1;
	
					}
					}
	
					else
					{
					data[ red_index] =  fluid.getDensity(x, y) * 255 / 6;
					data[ green_index] =  fluid.getDensity(x, y) * 215 / 6;
					}
				}
			}

}

//rainbow color
function vibgyor(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					//148, 0, 211 v
					else if(frames>=0 && frames<200)
					{
					data[ red_index] =  fluid.getDensity(x, y) * 148/ 6;
					data[ blue_index] =  fluid.getDensity(x, y) * 211 / 6;
					}
					//75, 0, 130 i
					else if (frames>=200 && frames<400)
					{
					data[ red_index] =  fluid.getDensity(x, y) * 75/ 6;
					data[ blue_index] =  fluid.getDensity(x, y) * 130 / 6;
					}
					else if (frames>=400 && frames<600)
					{
					data[ blue_index] =  fluid.getDensity(x, y) * 255/ 6;
					}
					else if (frames>=600 && frames<800)
					{
					
					data[ green_index] =  fluid.getDensity(x, y) * 255 / 6;
					}
					else if (frames>=800 && frames<1000)
					{
					data[ red_index] =  fluid.getDensity(x, y) * 255 / 6;
					data[ green_index] =  fluid.getDensity(x, y) * 255 / 6;
					}
					else if (frames>=1000 && frames<1200)
					{
					data[ red_index] =  fluid.getDensity(x, y) * 255 / 6;
					data[ green_index] =  fluid.getDensity(x, y) * 127 / 6;
					}
					else if (frames>=1200 && frames<=1400)
					{
					data[ red_index] =  fluid.getDensity(x, y) * 255 / 6;
					}
					else
					{
					frames=1;
					document.getElementById("reset").click();
					
					
					}
				}
			}

}

//blue color
function blue(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					
					data[blue_index] =  fluid.getDensity(x, y) * 255 / 6;
					
				}
			}

}


//pink color
function pink(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					
					
					
					data[red_index] =  fluid.getDensity(x, y) * 255 / 6;
					data[green_index] =  fluid.getDensity(x, y) * 20 / 6;
					data[blue_index] =  fluid.getDensity(x, y) * 147 / 6;
					
				}
			}

}


//white color

function white(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					
					data[blue_index] =  fluid.getDensity(x, y) * 255 / 6;
					data[green_index] =  fluid.getDensity(x, y) * 255 / 6;
					data[red_index] =  fluid.getDensity(x, y) * 255 / 6;
					
				}
			}

} 

//135,206,250 skyblue

function skyblue(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					data[red_index] =  fluid.getDensity(x, y) * 135 / 6;
					data[blue_index] =  fluid.getDensity(x, y) * 250 / 6;
					data[green_index] =  fluid.getDensity(x, y) * 206 / 6;
					
					
				}
			}

}

//	124-252-0 green

function green(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					
					data[red_index] =  fluid.getDensity(x, y) * 124 / 6;
					data[green_index] =  fluid.getDensity(x, y) * 252/ 6;
					
					
				}
			}

}


//	160-32-240 purple
function purple(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					
					data[red_index] =  fluid.getDensity(x, y) * 106 / 6;
					data[green_index] =  fluid.getDensity(x, y) * 32/ 6;
					data[blue_index] =  fluid.getDensity(x, y) * 240/ 6;
					
					
				}
			}

}


//red
function red(fluid,data,actual_frames,frames)
{
  
					
	for (var x = 0; x < width; x++) 
			{
				for (var y = 0; y < height; y++)
				{             
					var blue_index=4*(y * height + x)+2 ;
					var green_index=4*(y * height + x)+1 ;
	
					var red_index=4*(y * height + x) ;
					if(actual_frames>2500)
					{
   
					actual_frames=0;
					frames=1;
					document.getElementById("reset").click();
	 
	 
					}
					
					data[red_index] =  fluid.getDensity(x, y) * 255 / 6;
					
					
					
				}
			}

}


// draw density function

function draw_density(fluid)
{
	makeSecondCanvas(fluid);
	var context = canvas.getContext("2d");
	var width = secondCanvas.width;
	var height =secondCanvas.height;

	if (secondCanvasValues)
	 {
		var data = secondCanvasValues.data;
	   
		
		 {
			
					if(button_id=="cardinal")
					{
					
					cardinal(fluid,data,actual_frames,frames);
					}
					if(button_id=="blue")
					{
					
					blue(fluid,data,actual_frames,frames);
					}
					if(button_id=="skyblue")
					{
					
					skyblue(fluid,data,actual_frames,frames);
					}
					if(button_id=="green")
					{
					
					green(fluid,data,actual_frames,frames);
					}
					if(button_id=="purple")
					{
					
					purple(fluid,data,actual_frames,frames);
					}
					if(button_id=="red")
					{
					
					red(fluid,data,actual_frames,frames);
					}
					if(button_id=="vibgyor")
					{
					
					vibgyor(fluid,data,actual_frames,frames);
					}
					if(button_id=="white")
					{
					
					white(fluid,data,actual_frames,frames);
					}
					if(button_id=="pink")
					{
					
					pink(fluid,data,actual_frames,frames);
					}
				
		}
		context.putImageData(secondCanvasValues, 0, 0);
	} 
	else
	 {
		for (var x = 0; x < width; x++) 
		{
			for (var y = 0; y < height; y++) 
			{
				var d = fluid.getDensity(x, y) / 6;
				//context.fillStyle = 'red';
				context.fillStyle(0, 0, d, 1);
				context.fillRect(x, y, 1, 1);
			}
		}
	}
}


//velocity display function
function displayVelocity(fluid)
{
	var context = canvas.getContext("2d");
	context.save();
	context.lineWidth=1;
	var widthratio = canvas.width / secondCanvas.width;
	var heightratio = canvas.height / secondCanvas.height;
	
	
	context.fillStyle="white";
	context.fillRect(0, 0, canvas.width, canvas.height);
	context.strokeStyle = "rgb(0,0,255)";
	var vectorScale = 1;
	context.beginPath();
	for (var x = 0; x < secondCanvas.width; x++) 
	{
		for (var y = 0; y < secondCanvas.height; y++) 
		{
			context.moveTo(x * widthratio + 0.5 * widthratio, y * heightratio + 0.5 * heightratio);
			context.lineTo((x + 0.5 + vectorScale * fluid.getXVelocity(x, y)) * widthratio, 
						   (y + 0.5 + vectorScale * fluid.getYVelocity(x, y)) * heightratio);
		}
	}
	context.stroke();
	context.restore();
}

changeVisual = function(canvas) 
{
	
	canvas.width = 128;
	canvas.height = 128;
	return draw_density;
}
})();