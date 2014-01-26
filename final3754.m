                                          
function final3754()
    
    % Set domain constants

    xleft = 0.0;
    xright = 1.0;
    ybottom = 0.0;
    ytop = 1.0;
    
    h = 0.1;
    k = 0.1;
    
    nx = ((xright-xleft)/h) + 1;
    ny = ((ytop-ybottom)/k) + 1;

    n = nx*ny;
    
    nboundary = 2*(nx-1) + 2*(ny-1);
    
    
    % Set number of random walks for each point
    
    density = 100;  % (number of walks per boundary point)
    nwalks = density*nboundary;
    
    
    % Store boundary values in u
    
    u = zeros(n,1);
    
    for bndpt=2:nx
        xboundary = (bndpt-1)*h;
        yboundary = 0;
        u(bndpt,1) = 0;
    end
            
    for bndpt=2*nx:nx:nx*ny
        xboundary = (nx-1)*h;
        yboundary = ((bndpt/nx)-1)*k;
        u(bndpt,1) = 0.5 + 0.5* tanh(7.0*(yboundary - 0.5));
    end
            
    for bndpt=nx*ny-nx+1:nx*ny-1
        xboundary = (mod(bndpt,nx)-1)*h;
        yboundary = (ny-1)*k;
        u(bndpt,1) = 1.0;
    end
            
    for bndpt=1:nx:nx*ny-nx
        xboundary = 0.0;
        yboundary = floor(bndpt/nx)*k;
        u(bndpt,1) =0.5 + 0.5*(2*yboundary-1).^3;
    end
            
    
    % Set up x and y vectors (for plotting)
    
    for i=1:nx
        for j=1:ny
            point = nx*(j-1) + i;
            xpts(point,1) = (i-1)*h;
            ypts(point,1) = (j-1)*k;
        end
    end
    

    % Approximate solution for each interior point in the domain
    
    for istart=2:nx-1
        
        for jstart=2:ny-1
            
            % Set counts to 0 for all boundary points
            
            for pt=1:n
                count(pt,1) = 0;
            end
            
            % Perform n random walks, starting from (istart,jstart)
            
            for walk=1:nwalks  
         
                % Set x and y to starting point
                
                i = istart;
                j = jstart;
                
                % Walk through until you reach a boundary
                
                boundary = 0;
                
                while boundary==0
                    
                    % Set up probabilities for down, left, right, and up moves
                    % and set counter for point z to be 0

                    x = (i-1)*h;
                    y = (j-1)*k;
                    
                    probDown = 0.5*f(x,y)/(f(x,y) + (k*k)/(h*h));
                    probLeft = 0.5/(f(x,y)*((h*h)/(k*k)) + 1.0);
                    probRight = 0.5/(f(x,y)*((h*h)/(k*k)) + 1.0);
                    probUp = 0.5*f(x,y)/(f(x,y) + (k*k)/(h*h));
        
                    cumProb(1,1) = probDown;
                    cumProb(2,1) = cumProb(1,1) + probLeft;
                    cumProb(3,1) = cumProb(2,1) + probRight;
                    cumProb(4,1) = 1.0;
            
                    % Generate random number between 0 and 1
    
                    r = rand(1,1);
            
                    % Make move based on randon number
            
                    if r<cumProb(1,1)
                        j = j - 1;
                    elseif r<cumProb(2,1)
                        i = i - 1;
                    elseif r<cumProb(3,1)
                        i = i + 1;
                    else
                        j = j + 1;
                    end
             
                    % Check if at boundary, if so, add to count and break
                
                    if i==1 | i==nx | j==1 | j==ny
                        z = nx*(j-1) + i;
                        count(z,1) = count(z,1) + 1;
                        boundary = 1;
                    end
                    
                end
                
            end
            
            z = nx*(jstart-1) + istart;
            
            u(z,1) = 0.0;
            
            for bndpt=2:nx
                xboundary = (bndpt-1)*h;
                yboundary = 0;
                u(z,1) = u(z,1);
            end
            
            for bndpt=2*nx:nx:nx*ny
                xboundary = (nx-1)*h;
                yboundary = ((bndpt/nx)-1)*k;
                u(z,1) = u(z,1) + ((0.5 + 0.5* tanh(7.0*(yboundary - 0.5))))*(count(bndpt,1)/nwalks);
            end
            
            for bndpt=nx*ny-nx+1:nx*ny-1
                xboundary = (mod(bndpt,nx)-1)*h;
                yboundary = (ny-1)*k;
                u(z,1) = u(z,1) + 1*(count(bndpt,1)/nwalks);
            end
            
            for bndpt=1:nx:nx*ny-nx
                xboundary = 0.0;
                yboundary = floor(bndpt/nx)*k;
                u(z,1) = u(z,1) + (0.5 + 0.5*(2*yboundary-1).^3)*(count(bndpt,1)/nwalks);
            end
            
        end
        
    end
    
    
    % Make surface plot of approximate solution
    
    xplot = reshape(xpts,nx,ny);
    yplot = reshape(ypts,nx,ny);
    uplot = reshape(u,nx,ny)
    load ui.txt
    %u2 = u(:,any(u));
    error = norm(ui-u);     
    %surf(xplot,yplot,uplot,error)
    plot (error,u);

    
    
function z = f(x,y)         
    z = 1.0;
    
    
    
function z = g(x,y) 
        
    if abs(y-1.0)<0.001
        z = 1.0;
    else
        z = 0.0;
    end
     

