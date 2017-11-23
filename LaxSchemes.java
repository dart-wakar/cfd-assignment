import java.io.*;
import java.util.*;

public class LaxSchemes {

    double xmin,xmax,tmin,tmax,a = 1.0,dx,dt,sigma,x[],u0[],lfs[][],lws[][],exact[][];
    int num_x_nodes, num_t_nodes;

    public LaxSchemes(double xmin, double xmax, double tmax, double tmin, int num_x_nodes, int num_t_nodes){
        this.xmin = xmin;
        this.xmax = xmax;
        this.tmax = tmax;
        this.tmin = tmin;
        this.num_x_nodes = num_x_nodes;
        this.num_t_nodes = num_t_nodes;
        this.dx = (this.xmax - this.xmin) / (this.num_x_nodes - 1);
        this.dt = (this.tmax - this.tmin) / (this.num_t_nodes - 1);
        this.sigma = this.a * this.dt / this.dx;
    }

    public void generateMesh(){
        this.x = new double[this.num_x_nodes];
        this.u0 = new double[this.num_x_nodes];
        for(int i = 0; i < this.num_x_nodes; i++){
            this.x[i] = this.xmin + i*this.dx;
            this.u0[i] = Math.sin(2.0*Math.PI*this.x[i]);
        }
    }

    public void applyBoundaryConditions(){
        this.lfs = new double[this.num_t_nodes][this.num_x_nodes];
        this.lws = new double[this.num_t_nodes][this.num_x_nodes];
        for(int i = 0; i < this.num_x_nodes; i++){
            this.lfs[0][i] = Math.sin(2.0*Math.PI*this.x[i]);
            this.lws[0][i] = Math.sin(2.0*Math.PI*this.x[i]);
        }
        for(int i = 0; i < this.num_t_nodes; i++){
            this.lfs[i][0] = Math.sin((-2.0)*Math.PI*this.a*i*this.dt);
            this.lws[i][0] = Math.sin((-2.0)*Math.PI*this.a*i*this.dt);
            this.lfs[i][this.num_x_nodes-1] = Math.sin((-2.0)*Math.PI*this.a*i*this.dt);
            this.lws[i][this.num_x_nodes-1] = Math.sin((-2.0)*Math.PI*this.a*i*this.dt);
        }
    }

    public void solve(){
        this.exact = new double[this.num_t_nodes][this.num_x_nodes];
        for(int i = 1; i < this.num_t_nodes; i++){
            for(int j = 1; j < this.num_x_nodes - 1; j++){
                this.lfs[i][j] = ((this.lfs[i-1][j+1]+this.lfs[i-1][j-1])/2)-(this.sigma/2)*(this.lfs[i-1][j+1]-this.lfs[i-1][j-1]); // lax-friedrich
                this.lws[i][j]=this.lws[i-1][j]-((this.sigma/2)*(this.lws[i-1][j+1]-this.lws[i-1][j-1]))+((this.sigma*this.sigma/2)*(this.lws[i-1][j+1]-(2*this.lws[i-1][j])+this.lws[i-1][j-1])); // lax-wendroff
                this.exact[i][j] = Math.sin(2*Math.PI*((this.xmin+j*this.dx)-(this.a*(0+i*this.dt))));
            }
        }
    }

    public void printLastRow(){
        for(int i = 0; i < this.num_x_nodes; i++){
            System.out.println(this.lfs[this.num_t_nodes - 1][i]);
            System.out.println(this.lfs[this.num_t_nodes - 1][i]);
        }
    }

    public void writeIntoFile(){
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter("lax_friedrich.txt"));
            for(int i  = 0; i < this.num_t_nodes; i++){
                for(int j = 0; j < this.num_x_nodes; j++){
                    bw.write(this.lfs[i][j] + ((j == this.num_x_nodes - 1) ? "" : " "));
                }
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e){
            System.out.println("Error in writing into file " + e);
        }
        try {
            BufferedWriter bw1 = new BufferedWriter(new FileWriter("lax_wendroff.txt"));
            for(int i  = 0; i < this.num_t_nodes; i++){
                for(int j = 0; j < this.num_x_nodes; j++){
                    bw1.write(this.lws[i][j] + ((j == this.num_x_nodes - 1) ? "" : " "));
                }
                bw1.newLine();
            }
            bw1.flush();
        } catch (IOException e){
            System.out.println("Error in writing into file " + e);
        }
        try {
            BufferedWriter bw2 = new BufferedWriter(new FileWriter("exact.txt"));
            for(int i  = 0; i < this.num_t_nodes; i++){
                for(int j = 0; j < this.num_x_nodes; j++){
                    bw2.write(this.exact[i][j] + ((j == this.num_x_nodes - 1) ? "" : " "));
                }
                bw2.newLine();
            }
            bw2.flush();
        } catch (IOException e){
            System.out.println("Error in writing into file " + e);
        }
    }

    public static void main(String args[]){
        LaxSchemes laxSchemes = new LaxSchemes(0.0, 1.0, 1.0, 0.0, 100, 100);
        laxSchemes.generateMesh();
        laxSchemes.applyBoundaryConditions();
        laxSchemes.solve();
        laxSchemes.printLastRow();
        laxSchemes.writeIntoFile();
    }

}