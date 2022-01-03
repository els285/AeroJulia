import json
from os import stat
# from os import XATTR_SIZE_MAX 
import plotly.graph_objects as go
import numpy as np
from dataclasses import dataclass




class Shock_Expansion_Plot:

    """
    Shock_Expansion_Plot generates a plotly figure where the supersonic body and various flow regions are drawn
    """

    def __init__(self,filename):

        # Initialise figure
        self.fig = go.Figure()

        # Extract data from JSON file
        with open(filename) as f:
            self.data = json.load(f)

        # Define paramters from JSON data
        self.body_data = self.data['Body']["profiles"]
        self.FlowSolution = self.data["FlowSolution"]

        # Compute extrema
        self.extrema = self.maximal_extent()


    def draw_plot(self):

        self.plot_body()     
        self.draw_upstream("upper")
        self.draw_upstream("lower")
        self.draw_other_regions("upper")
        self.draw_other_regions("lower")
        self.draw_downstream("upper")
        self.draw_downstream("lower")
        return self.fig



    def maximal_extent(self):

        """
        The extrema of the plot are necessary for drawing the shock region polygons, which are defined by shock or exapnsion fan angles
        """

        @dataclass
        class Extrema:
            Xmin: float
            Xmax: float 
            Ymin: float 
            Ymax: float

        x_points = set([p['x'] for p in self.body_data['upper']] + [p['x'] for p in self.body_data['lower']])
        y_points = set([p['y'] for p in self.body_data['upper']] + [p['y'] for p in self.body_data['lower']])

        x_min , x_max  =  min(x_points) , max(x_points)
        y_min , y_max  =  min(y_points) , max(y_points) 

        # Plot should be square and twice the length of the largest dimension
        y_scale = 3
        x_scale_upstream = 10
        x_scale_downstream = 20
        maxX = max(abs(x_max),abs(x_min))
        maxY = max(abs(y_max),abs(y_min))
        plot_Xmin , plot_Xmax  = x_min - x_scale_upstream*maxY , x_max + x_scale_downstream*maxY

        plot_Ymin , plot_Ymax  = y_min - y_scale*maxY , y_max + y_scale*maxY

        return Extrema(plot_Xmin,plot_Xmax,plot_Ymin,plot_Ymax)



    def plot_body(self):

        """ Plots the supersonic body """

        self.fig.add_trace(go.Scatter(x=[p['x'] for p in self.body_data['upper']],
                                y=[p['y'] for p in self.body_data['upper']],  
                                fillcolor="RoyalBlue", fill="toself"    , line_width=0))

        self.fig.add_trace(go.Scatter(x=[p['x'] for p in self.body_data['lower']],
                                y=[p['y'] for p in self.body_data['lower']],  
                                fillcolor="RoyalBlue", fill="toself"    , line_width=0))

    @staticmethod
    def get_x_from_y(angle,point0,y1):

        """ 
        Computing the final parameter of a straigh line
        Takes the initial point, which is on the body surface 
        """

        x0 = point0[0]
        y0 = point0[1]
        y1 = y1
        return (y1 - y0)/np.tan(angle) + x0


    def draw_upstream(self,edge): 

        """
        For drawing the initial flow state (in two separate parts) before any body-interaction
        """

        edge_profile  =  self.body_data[edge]
        f = self.FlowSolution[edge][1]

        if "incoming_angle" in f.keys() and "outgoing_angle" in f.keys():
            angle = f["incoming_angle"]
        elif "Cone_Angle" in f.keys() and "Shock_Angle" in f.keys():
            angle = f["Shock_Angle"]
        # Make upstream polygon
        if edge=="upper":
            max_point = self.extrema.Ymax
        elif edge=="lower":
            max_point = self.extrema.Ymin
            angle = -angle

        # Define points of quadrilateral
        x0,y0 = edge_profile[0]["x"] , edge_profile[0]["y"]
        x1,y1 = self.get_x_from_y(angle,(x0,y0),max_point),max_point
        x2,y2 = self.extrema.Xmin,max_point
        x3,y3 = self.extrema.Xmin,edge_profile[0]["y"]

        # Draw
        self.fig.add_trace(go.Scatter(x=[x0,x1,x2,x3,x0] , y=[y0,y1,y2,y3,y0],fillcolor="Red",fill="toself"    , line_width=0))




    def draw_Prantl_Meyer(self,p0,p1,incoming_angle,outgoing_angle,next_angle,Downstream_flow,edge):

        """
        Draws both the expansion fan and the 
        """

        if edge=="upper":
            max_point = self.extrema.Ymax
        elif edge=="lower":
            max_point = self.extrema.Ymin
            incoming_angle = -incoming_angle
            outgoing_angle = -outgoing_angle
            next_angle     = -next_angle

        # Draw the fan
        # Proceed anti-clockwise
        x0,y0 = p0["x"],p0["y"]
        x1,y1 = self.get_x_from_y(outgoing_angle,(x0,y0),max_point) , max_point
        x2,y2 = self.get_x_from_y(incoming_angle,(x0,y0),max_point) , max_point

        self.fig.add_trace(go.Scatter(x=[x0,x1,x2,x0] , y=[y0,y1,y2,y0],fillcolor="Orange",fill="toself"    , line_width=0))


        # Draw the downstream quadrilateral
        xA,yA = p1["x"],p1["y"]
        xB,yB = self.get_x_from_y(next_angle,(xA,yA),max_point) , max_point
        self.fig.add_trace(go.Scatter(x=[x0,xA,xB,x1,x0] , y=[y0,yA,yB,y1,y0],fillcolor="Green",fill="toself"    , line_width=0))


    def draw_oblique_shock(self,p0,p1,this_shock_angle,next_angle,Downstream_flow,edge):

        if edge=="upper":
            max_point = self.extrema.Ymax
        elif edge=="lower":
            max_point = self.extrema.Ymin
            this_shock_angle = -this_shock_angle
            next_angle       = -next_angle


        x0,y0 = p0["x"],p0["y"]
        x1,y1 = p1["x"],p1["y"]
        x2,y2 = self.get_x_from_y(next_angle,(x1,y1),max_point) , max_point
        x3,y3 = self.get_x_from_y(this_shock_angle,(x0,y0),max_point) , max_point

        self.fig.add_trace(go.Scatter(x=[x0,x1,x2,x3,x0] , y=[y0,y1,y2,y3,y0],fillcolor="Purple",fill="toself"    , line_width=0,hovertext="Jgjhhjg"))



    def draw_other_regions(self,edge):

        edge_profile  =  self.body_data[edge]
        flow_solution =  self.FlowSolution[edge]        

        for i,p in enumerate(edge_profile[:-1]):
            p0 = edge_profile[i]
            p1 = edge_profile[i+1]
            f = flow_solution[i+1]
            f_next = flow_solution[i+2]
            
            # Set next angle for 
            if "incoming_angle" in f_next.keys():    next_angle = f_next["incoming_angle"]
            elif "Shock_Angle" in f_next.keys():     next_angle = f_next["Shock_Angle"] 

            # Draw flow polygon
            if "incoming_angle" in f.keys() and "outgoing_angle" in f.keys():
                self.draw_Prantl_Meyer(p0,p1,f["incoming_angle"],f["outgoing_angle"],next_angle,f["Downstream"],edge)

            elif "Cone_Angle" in f.keys() and "Shock_Angle" in f.keys():
                self.draw_oblique_shock(p0,p1,f["Shock_Angle"],next_angle,f["Downstream"],edge)



    def draw_downstream(self,edge):

        edge_profile  =  self.body_data[edge]
        flow_solution =  self.FlowSolution[edge]    

        p = edge_profile[-1]
        f = flow_solution[-1]


        if "incoming_angle" in f.keys() and "outgoing_angle" in f.keys():
            angle = f["incoming_angle"]
        elif "Cone_Angle" in f.keys() and "Shock_Angle" in f.keys():
            angle = f["Shock_Angle"]

        if edge=="upper":
            max_point = self.extrema.Ymax
        elif edge=="lower":
            max_point = self.extrema.Ymin
            angle = -angle

        x0,y0 = p["x"] , p["y"]
        x1,y1 = self.extrema.Xmax , p["y"]
        x2,y2 = self.extrema.Xmax , max_point
        x3,y3 = self.get_x_from_y(angle,(x0,y0),max_point),max_point

        self.fig.add_trace(go.Scatter(x=[x0,x1,x2,x3,x0] , y=[y0,y1,y2,y3,y0],fillcolor="Red",fill="toself"    , line_width=0))




def main():
    plot=Shock_Expansion_Plot("foo.json")
    fig = plot.draw_plot()
    fig.show()

main()