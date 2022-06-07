#include "Format.h"
#include "spectra_v2_ampt.h"

using namespace std;

vector<double> Decay(vector<double> &pin, double m00);

void spectra_v2_ampt(int syst, int pr, int centr) {

    bool if_check_Nch = false;
    TH1D *hist_Nch = new TH1D("hist_Nch","hist_Nch",500,-0.5,499.5);

    bool pions = false;
    bool real_spectra_phi = true;
    double eta_cut = 0.5;
    bool plot_v2 = false;

    ///****************inizializing*******************************
    gStyle->SetLineWidth(3);

    double px,py,pz,m,y,pt,p_tot, nch, v2, phi_bbc,phi_cnt, phi_fvtx, ncoll_ampt_read, bimp;
    double nlines = 0,nlines_pi0 = 0, ncollisions = 0, ncollisions_pi0 = 0, n_bad =0, ncoll_ampt=0;
    TFile *f;
    TH1D *pt_dist_SS = new TH1D("pt_dist_SS", "pT distribution SS", 12, pt_bins);
    TH1D *y_dist_SS = new TH1D("y_dist_SS", "Y distribution SS", 150, -3, 3);
    TH1D *res[3];
    for (int i = 0; i < 3; ++i) {
        sprintf(name,"res_%d",i);
        res[i] = new TH1D(name, name, 200, -1,1);
    }
    TH2D *v2_SS = new TH2D("v2_dist_SS", "#varphi v_{2} p+Au", 2000, -1, 1, 50, 0, 5);
    TH1D *pt_dist_SS_pi0 = new TH1D("pt_dist_SS_pi0", "pT distribution SS_pi0", 12, pt_bins);
    TNtuple *ntuple, *ntuple_pi0;
    double ratio[12], ratio_er[12], pt_r[12], pt_r_pi0[12], dpt_r_pi0[12];

    TString dir = gSystem->UnixPathName(__FILE__);
    ifstream in, in1;
    char names[10][200];

    double Ncollis = Ncoll[syst/2][centr];
    if(pr==0)
    {
        sprintf(names[0],"output_code/phi_%s_%s.root",syst_names[syst/2],ampt_type[syst%2]);
        f = new TFile(names[0], "RECREATE");
        ntuple = new TNtuple("ntuple", "ntuple", "nch:px:py:pz:phi");
        TH1D *h1 = new TH1D("h1", "x distribution", 5000, 0, 5000);
        TH2D *h2 = new TH2D("h2", "x distribution", 5000, 0, 5000,500,0,500);
        char filename[200],  filenamepp[200];
        for (int i = 0; i < 1000; ++i) {
            //sprintf(filename,"output_taxi/new/SS_30mb_35_nosh_q/output_taxi/ampt_%d.dat",i);
           // sprintf(filename,"/home/yoren/Documents/AMPT/backup/SS_sm_30mb_3040/output_taxi/ampt_%d.dat",i);
            sprintf(filename,"/home/yoren/bnl/AMPT/work/new/%s_30mb_3050_%s/output_taxi/ampt_%d.dat",syst_names[syst/2],ampt_type[syst%2],i);
            if(if_check_Nch)sprintf(filename,"/home/yoren/bnl/AMPT/work/new_test/pAl_sm_test/output_taxi/ampt_%d.dat",i);

            in.open(filename);
            while (1) {
                in >> px;
                if (!in.good()) break;
                if(px>=100)
                {
                    in >>px>>py>>bimp>>pz>>m>>px>>px>>ncoll_ampt_read>>px>>phi_bbc>>phi_cnt>>phi_fvtx;
                    if (nlines < 10) printf("Nch=%f, phi = %f, %f, %f\n", py,phi_bbc,phi_cnt,phi_fvtx);
                    if(py<0) continue;

                    if(if_check_Nch) hist_Nch->Fill(px);
                    //cout<<py<<endl;
                    //if(fabs(ncoll_ampt_read-320)>10)continue;
                    nch = py;
                    h1->Fill(nch);
                    h2->Fill(nch,m+pz);
                    if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]])
                    {
                        ncoll_ampt+=ncoll_ampt_read;
                        res[0]->Fill(cos(2*(phi_fvtx-phi_cnt)));
                        res[1]->Fill(cos(2*(phi_bbc-phi_fvtx)));
                        res[2]->Fill(cos(2*(phi_bbc-phi_cnt)));
                        ncollisions++;
                    }
                }
                else
                {
                    in >> py >> pz >> m;
                    if (nlines < 50) printf("px=%f, py=%f, pz=%f, m=%f\n", px,py,pz,m);
                    pt=sqrt(px*px+py*py);
                    if(pt<0.8) continue;
                    p_tot=sqrt(pt*pt+pz*pz);//+m*m
                    y=0.5*log((p_tot+pz)/(p_tot-pz));
                    if(fabs(y)>=eta_cut) continue;
                    double phi = atan2(py,px);
                    v2 = cos(2*(phi-phi_bbc));
                    //v2= (px*px-py*py)/pt/pt;
                    if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) v2_SS->Fill(v2,pt);
                    if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) pt_dist_SS->Fill(pt);
                    if (pt>1.0&&pt<4.0) y_dist_SS->Fill(y);
                    ntuple->Fill(nch,px,py,pz,phi_bbc);
                    nlines++;
                }
            }
            in.close();
        }
        if(pions){
            TH1D *h11 = new TH1D("h11", "x distribution", 5000, 0, 5000);
            ntuple_pi0 = new TNtuple("ntuple_pi0", "ntuple_pi0", "nch:pt");
            for (int i = 0; i < 1000; ++i) {
                //sprintf(filename,"output_taxi/SS_pi0/ampt_%d.dat",i);
                sprintf(filename,"/home/yoren/Documents/AMPT/backup/%s_pi0/output_taxi/ampt_%d.dat",syst_names[syst/2],i);
                in1.open(filename);
                while (1) {
                    in1 >> px;
                    if (!in1.good()) break;
                    if(px>=100)
                    {
                        in1>>px>>py>>bimp>>pz>>m>>px>>px>>ncoll_ampt_read>>py>>phi_bbc>>phi_cnt>>phi_fvtx;
                        if(py<0) continue;
                        if (nlines_pi0 < 10) printf("Nch=%f\n", py);
                        h11->Fill(py);
                        nch = py;
                        if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) ncollisions_pi0++;
                    }
                    else
                    {
                        in1 >> py >> pz >> m;
                        if(fabs(m-0.137)>0.0031)
                        {
                            vector <double> pin = {px, py, pz};
                            vector <double> ppp = Decay(pin,m);
                            if (ppp[0]!=0&&ppp[1]!=0) {px = ppp[0]; py = ppp[1]; pz = ppp[2];}
                            //else cout<<m<<endl;
                            //cout<<sqrt(px*px+py*py)/sqrt(pin[0]*pin[0]+pin[1]*pin[1])<<endl;
                            //cout<<"WTF1:"<<m<<endl; n_bad++; break;
                        }
                        if (nlines_pi0 < 50) printf("px=%f, py=%f, pz=%f, m=%f\n", px,py,pz,m);
                        pt=sqrt(px*px+py*py);
                        if(pt<0.8) continue;
                        p_tot=sqrt(pt*pt+pz*pz);//+m*m
                        y=0.5*log((p_tot+pz)/(p_tot-pz));
                        if(fabs(y)>=eta_cut) continue;
                        if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) pt_dist_SS_pi0->Fill(pt);
                        ntuple_pi0->Fill(nch,pt);
                        nlines_pi0++;
                    }
                }
                in1.close();
            }

            cout<<ncollisions_pi0<<"  "<<ncollisions<<endl;
            /*pt_dist_SS->Divide(pt_dist_SS_pi0);
            pt_dist_SS->Scale(ncollisions_pi0/ncollisions);
            pt_dist_SS->Draw();*/

        }

        cout<<ncollisions<<"  "<<ncollisions_pi0<<endl;
        printf("found %0.f collisions and %0.f points\n", ncollisions_pi0,n_bad);
        f->Write();
        if(if_check_Nch) hist_Nch->Write();
    }

    if(pr>0)
    {
        sprintf(names[1],"output_code/phi_%s_%s.root",syst_names[syst/2],ampt_type[syst%2]);
        f = TFile::Open(names[1]);
        TH1D *h2 = (TH1D*) gROOT->FindObject("h1");
        TH2D *h3 = (TH2D*) gROOT->FindObject("h2");
        ntuple = (TNtuple*) gROOT->FindObject("ntuple");
        ncollisions = 0;
        for (int i = h2->FindBin(nch_bin[syst][ncentr_bin[syst/2][2*centr]]); i < h2->FindBin(nch_bin[syst][ncentr_bin[syst/2][2*centr+1]])+1; ++i) {
            ncollisions+=h2->GetBinContent(i);
        }
        TH1D *h4 = h3->ProjectionY("h4",(int) (nch_bin[syst][ncentr_bin[syst/2][2*centr]]+0.5),(int) (nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]+0.5));
        double npart = h4->GetMean();
        Ncollis =Ncoll[syst/2][centr];// 1.616*npart+0.004614*npart*npart;
        cout<<Ncollis<<" "<<Ncoll[centr]<<endl;
        //ncollisions= 1830;
        for (int ii = 0; ii < ntuple->GetEntries(); ++ii) {
            ntuple->GetEntry(ii);
            nlines++;
            float *time = ntuple->GetArgs();
            nch = (double) time[0];
            px = (double) time[1];
            py = (double) time[2];
            pz = (double) time[3];
            phi_bbc = (double) time[4];
            pt=sqrt(px*px+py*py);
            double phi = atan2(py,px);
            v2 = cos(2*(phi-phi_bbc));
            //v2 = (px*px-py*py)/pt/pt;
            if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) v2_SS->Fill(v2,pt);
            if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) pt_dist_SS->Fill(pt);
        }
        if(pions){
            ntuple_pi0 = (TNtuple*) gROOT->FindObject("ntuple_pi0");
            ncollisions_pi0 = 0;
            TH1D *h11 = (TH1D*) gROOT->FindObject("h11");
            for (int i = h11->FindBin(nch_bin[syst][ncentr_bin[syst/2][2*centr]]); i < h11->FindBin(nch_bin[syst][ncentr_bin[syst/2][2*centr+1]])+1; ++i) {
                ncollisions_pi0+=h11->GetBinContent(i);
            }
            for (int ii = 0; ii < ntuple_pi0->GetEntries(); ++ii) {
                ntuple_pi0->GetEntry(ii);
                nlines++;
                float *time = ntuple_pi0->GetArgs();
                nch = (double) time[0];
                pt = (double) time[1];
                if (nch>nch_bin[syst][ncentr_bin[syst/2][2*centr]]&&nch<nch_bin[syst][ncentr_bin[syst/2][2*centr+1]]) pt_dist_SS_pi0->Fill(pt);
            }

            cout<<ncollisions<<"  "<<ncollisions_pi0<<endl;
        }
    }

//pt_dist_SS->Draw();

    int Nevnt = ncollisions;//140000
    cout<<"Ncoll: "<<  ncoll_ampt/ncollisions<<" "<<Ncoll[centr]<<endl;

    if(!pions){

        double yeilds_pythia[N], yeilds_pythia_err[N], pT[N], dpT[N], R_SS[N], R_SS_er[N];

        for (int i = 0; i < N; i++) {
            yeilds_pythia[i]=pt_dist_SS->GetBinContent(i+1);
            yeilds_pythia_err[i] = sqrt(yeilds_pythia[i]);
            pT[i] = pt_dist_SS->GetBinCenter(i+1);
            dpT[i] = pt_dist_SS->GetBinWidth(i+1)/2;
            R_SS[i]=yeilds_pythia[i]/pT[i]/(2*dpT[i])/(2*3.14159)/Nevnt/yo_pp[i]/Ncollis/(2*eta_cut);
            R_SS_er[i]=R_SS[i]*sqrt(pow(yeilds_pythia_err[i]/yeilds_pythia[i], 2) + pow(eyo_pp[i]/yo_pp[i], 2));
            cout<<i<<" "<<2*dpT[i]<<" "<<pT[i]<<" "<<R_SS[i]<<" "<<Ncollis<<" "<<Nevnt<<" "<<yeilds_pythia[i]<<" "<<yo_pp[i]<<endl;
        }

        cout<<"pT["<<N-1<<"]={";
        for (int j = 0; j < N-1; ++j) {
            if(j!=9)cout<<pT[j]<<", ";
        }
        cout<<pT[N-1]<<"};"<<endl;
        cout<<"dpT["<<N-1<<"]={";
        for (int j = 0; j < N-1; ++j) {
            if(j!=9)cout<<dpT[j]<<", ";
        }
        cout<<dpT[N-1]<<"};"<<endl;
        cout<<"R_"<<syst_names[syst/2]<<"_"<<ampt_type[syst%2]<<"["<<N-1<<"]={";
        for (int j = 0; j < N-1; ++j) {
            if(j!=9)cout<<R_SS[j]<<", ";
        }
        cout<<R_SS[N-1]<<"};"<<endl;
        cout<<"R_stat_"<<syst_names[syst/2]<<"_"<<ampt_type[syst%2]<<"["<<N-1<<"]={";
        for (int j = 0; j < N-1; ++j) {
            if(j!=9)cout<<R_SS_er[j]<<", ";
        }
        cout<<R_SS_er[N-1]<<"};"<<endl;

        TGraphErrors *graph = new TGraphErrors(N,pT,R_SS,dpT,R_SS_er);
        TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
        c2->SetLeftMargin(0.13);
        c2->SetBottomMargin(0.13);
        c2->SetTickx();
        c2->SetTicky();
        char ord_name[200];
        sprintf(ord_name,"R_{%s}",syst_names[syst/2]);
        Format_Pad(0.5,8.2,0.2,2.2,"p_{T}, GeV/c",ord_name,1,1.05,0.055,0.05,"");
        Format_Graph(graph,22,3,2,2,2,1,1);
        //Format_Graph(graph,22,0,0,0,0,1,1);
        if(pr<2)graph->Draw("P");

        TGraphErrors *graph_real = new TGraphErrors(N,pt_phi,raa_SS[centr],pt_er1,rae_SS[centr]);
        Format_Graph(graph_real,22,3,4,2,4,1,1);
        TGraphErrors *graph_real_s = new TGraphErrors(N,pt_phi,raa_SS[centr],pt_er1,ras_SS[centr]);
        Format_Graph(graph_real_s,22,0,4,2,4,1,1);
        graph_real_s->SetFillColorAlpha(4,0.5);
        //if(pr>1) Format_Graph(graph,22,0,0,0,0,1,1);
        if(pr<2)graph_real->Draw("P");
        if(pr>1)graph_real_s->Draw("2P");



        if(pr>1)
        {
            graph->SetFillColorAlpha(2,0.3);
            graph->RemovePoint(9);
            graph->Draw("same e3");
            graph_real->Draw("sameP");
        }

        TLegend *legend = new TLegend(0.55,0.45,0.85,0.65);
        legend->SetFillColorAlpha(0,0);
        legend->SetLineWidth(2);
        legend->SetTextSize(0.06);
        legend->AddEntry(graph_real,"real data","lep");
        legend->AddEntry(graph,"AMPT","f");
        //legend->Draw();

        TLegend *legend2 = new TLegend(0.2,0.8,0.85,0.9);
        sprintf(name,"#varphi, %d-%d%% %s",central_bin[centr*2],central_bin[centr*2+1],syst_titles[syst/2]);
        legend2->SetHeader(name,"c");
        legend2->SetLineWidth(0);
        legend2->SetTextSize(0.06);
        legend2->SetFillColorAlpha(0,0);
        legend2->Draw();

        TLegend *legend3 = new TLegend(0.2,0.7,0.85,0.8);
        legend3->SetHeader("|y|<0.35, #sqrt{s_{NN}}=200 GeV","c");
        legend3->SetLineWidth(0);
        legend3->SetTextSize(0.06);
        legend3->SetFillColorAlpha(0,0);
        legend3->Draw();

       // DrawLogo(0.12,0,0.5,0.33);

        sprintf(name,"output_code/%s_%d%d.png",syst_names[syst/2], central_bin[centr*2],central_bin[centr*2+1]);
        c2->SaveAs(name);

        TProfile *v2_hist = v2_SS->ProfileY();

        if(plot_v2)
        {
            TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
            gStyle->SetOptStat(0000);

            c3->SetLeftMargin(0.13);
            c3->SetBottomMargin(0.13);
            c3->SetTickx();
            c3->SetTicky();
            sprintf(name,"#varphi#rightarrowK^{+}K^{-}, %d-%d%% %s, 200 GeV",central_bin[2*centr],central_bin[2*centr+1],syst_titles[syst]);
            Format_Pad(0.,5.,0.002,0.35,"p_{T}, GeV/c","V_{2}  ",1,1,0.055,0.045,name);
            v2_hist->Rebin(5);
            Format_Profile(v2_hist,21,3,2,2,2,1,1);
            if(pr>1) Format_Profile(v2_hist,21,0,2,0,2,1,1);
            if(pr<2)v2_hist->Draw("same");
            if(pr>1)
            {
                v2_hist->SetFillColorAlpha(2,1);
                v2_hist->SetFillStyle(3001);
                v2_hist->Draw("sameE4");
            }
            TLegend *legend1 = new TLegend(0.2,0.65,0.45,0.85);
            legend1->SetFillColorAlpha(0,0);
            legend1->SetLineWidth(2);
            legend1->AddEntry(graph_real,"real data","lep");
            legend1->AddEntry(v2_hist,"AMPT","f");
            legend1->Draw();
            sprintf(name,"output_code/v2_%s_AMPT_%d%d.png",syst_names[syst/2],central_bin[centr*2],central_bin[centr*2+1]);
            c3->SaveAs(name);
        }

        cout<<"{";
        for (int j = 1; j < v2_hist->GetNbinsX(); ++j) {
            cout<<v2_hist->GetBinContent(j)<<", ";
        }
        cout<<v2_hist->GetBinContent(v2_hist->GetNbinsX())<<"}"<<endl;
        cout<<"{";
        for (int j = 1; j < v2_hist->GetNbinsX(); ++j) {
            cout<<v2_hist->GetBinError(j)<<", ";
        }
        cout<<v2_hist->GetBinError(v2_hist->GetNbinsX())<<"}"<<endl;
        cout<<"{";
        for (int j = 1; j < v2_hist->GetNbinsX(); ++j) {
            cout<<v2_hist->GetBinCenter(j)<<", ";
        }
        cout<<v2_hist->GetBinCenter(v2_hist->GetNbinsX())<<"}"<<endl;

    } else{
        for (int ibin = 0; ibin < n_pt_sp; ++ibin) {
            pt_r_pi0[ibin]=pt_dist_SS_pi0->GetBinCenter(pi0_bins[ibin]);
            dpt_r_pi0[ibin]=pt_dist_SS_pi0->GetBinWidth(pi0_bins[ibin]);
        }

        for (int ibin = 0; ibin < n_pt_sp; ++ibin) {
            ratio[ibin] = spectra_SS[centr][ibin]/pt_dist_SS_pi0->GetBinContent(pi0_bins[ibin])*pt_r_pi0[ibin]*dpt_r_pi0[ibin]* ncollisions_pi0*2*3.14159*(2*eta_cut);
            ratio_er[ibin]=ratio[ibin]*sqrt(spectra_SS_stat[centr][ibin]+1./pt_dist_SS_pi0->GetBinContent(pi0_bins[ibin]));
            cout<<pt_r_pi0[ibin]<<" "<<spectra_SS[centr][ibin]<<" "<<pt_dist_SS_pi0->GetBinContent(pi0_bins[ibin])/pt_r_pi0[ibin]/dpt_r_pi0[ibin]/ncollisions_pi0/2/3.14159<<endl;
            dpt_r_pi0[ibin]/=2;
        }
        if(!real_spectra_phi){
            for (int ibin = 0; ibin < pt_dist_SS->GetNbinsX(); ++ibin) {
                pt_r_pi0[ibin]=pt_dist_SS_pi0->GetBinCenter(ibin+1);
                dpt_r_pi0[ibin]=pt_dist_SS_pi0->GetBinWidth(ibin+1)/2;
                ratio[ibin] = pt_dist_SS->GetBinContent(ibin+1)/pt_dist_SS_pi0->GetBinContent(ibin+1)/ncollisions* ncollisions_pi0;
                ratio_er[ibin]=ratio[ibin]*sqrt(1./pt_dist_SS->GetBinContent(ibin+1)+.1/pt_dist_SS_pi0->GetBinContent(ibin+1));
            }
        }

        for (int i = 0; i < n_pt_sp; ++i) {
            cout<<ratio[i]<<", ";
        }
        cout<<endl;
        for (int i = 0; i < n_pt_sp; ++i) {
            cout<<ratio_er[i]<<", ";
        }
        cout<<endl;

        TGraphErrors *graph = new TGraphErrors(n_pt_sp,pt_r_pi0,ratio,dpt_r_pi0,ratio_er);
        if(!real_spectra_phi)TGraphErrors *graph = new TGraphErrors(12,pt_r_pi0,ratio,dpt_r_pi0,ratio_er);
        TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
        c2->SetLeftMargin(0.13);
        c2->SetBottomMargin(0.13);
        c2->SetTickx();
        c2->SetTicky();
        char ord_name1[200];
        sprintf(ord_name1,"R_{%s}",syst_names[syst/2]);
        Format_Pad(0.5,8.2,0.0,0.6,"p_{T}, GeV/c",ord_name1,1,1.05,0.055,0.05,"");
        Format_Graph(graph,22,3,2,2,2,1,1);
        //Format_Graph(graph,22,0,0,0,0,1,1);
        if(pr<2)graph->Draw("P");
    }

    cout<<sqrt(res[1]->GetMean()*res[2]->GetMean()/res[0]->GetMean())<<endl;

    if(if_check_Nch)cout<<hist_Nch->GetMean()<<"  "<<hist_Nch->GetEntries()<<endl;


}


vector<double> Decay(vector<double> &pin, double m00){

    vector<double> pp_rot ={0,0,0,0,0,0};

    double m01 = 0, m02 =0;

    double rand = f1->GetRandom();
    if((m00==0.494||m00==0.498)&&rand<=0.614){ m01 = 0.135, m02 =0.135;}//kstort klong
    else if((m00==1.116)&&rand<=0.36*1.2){ m01 = 0.135, m02 =0.938;}//lambda0
    else if((m00==1.189)&&rand<=0.52){ m01 = 0.135, m02 =0.937;}//Sigma+
    else if((m00==1.192)&&rand<=1){ cout<<"yolo"<<endl;}//Sigma0
    else if((m00==1.315)&&rand<=1.){ m01 = 0.135, m02 =1.116;}//Xi- lamda else
    else if((m00==2.453)&&rand<=1.){ m01 = 0.135, m02 =2.286;}//Xic- lamda else
    else if((m00==2.452)&&rand<=1.){ m01 = 0.135, m02 =2.286;}//Xic- lamda else
    else if((m00==1.672)&&rand<=0.08){ m01 = 0.135, m02 =1.189;}//Xic- lamda else 0.08
        //else if(m00==1.197) { m01 = 0.135, m02 =0.135;}//bred
        //else if(m00==1.321) { m01 = 0.135, m02 =0.135;}//bred
        //else if(m00==2.452) { m01 = 0.135, m02 =0.135;}//bred
        //else if(m00==2.285) { m01 = 0.135, m02 =0.135;}//bred
    else return pp_rot;



    double p11 = sqrt(pow(m01,4)-2*pow(m01*m02,2)-2*pow(m01*m00,2)+pow(m02,4)-2*pow(m02*m00,2)+pow(m00,4))/(2*m00);

    double theta = f2->GetRandom();
    double phi = f3->GetRandom();

    double p_dec[3] = {p11*cos(phi)*sin(theta),p11*sin(phi)*sin(theta),p11*cos(theta)};

    double ptot = sqrt(pin[0]*pin[0]+pin[1]*pin[1]+pin[2]*pin[2]);
    double E = sqrt(m00*m00+ptot*ptot);
    double beta = ptot/E;
    double gamma = E/m00;

    double theta_in = atan2(sqrt(pin[0]*pin[0]+pin[1]*pin[1]),pin[2]);
    double phi_in = atan2(pin[1],pin[0]);

    double E01 = sqrt(m01*m01+p11*p11);
    double E02 = sqrt(m02*m02+p11*p11);

    vector<double> pp_daugters;

    pp_daugters.push_back(p_dec[0]);
    pp_daugters.push_back(p_dec[1]);
    pp_daugters.push_back((p_dec[2]+beta*E01)*gamma);
    pp_daugters.push_back(-p_dec[0]);
    pp_daugters.push_back(-p_dec[1]);
    pp_daugters.push_back((-p_dec[2]+beta*E02)*gamma);


    pp_rot[0]=pp_daugters[0]*cos(theta_in)*cos(phi_in)+pp_daugters[2]*sin(theta_in)*cos(phi_in)-pp_daugters[1]*sin(phi_in);
    pp_rot[3]=pp_daugters[3]*cos(theta_in)*cos(phi_in)+pp_daugters[5]*sin(theta_in)*cos(phi_in)-pp_daugters[4]*sin(phi_in);
    pp_rot[1]=pp_daugters[0]*cos(theta_in)*sin(phi_in)+pp_daugters[2]*sin(theta_in)*sin(phi_in)+pp_daugters[1]*cos(phi_in);
    pp_rot[4]=pp_daugters[3]*cos(theta_in)*sin(phi_in)+pp_daugters[5]*sin(theta_in)*sin(phi_in)+pp_daugters[4]*cos(phi_in);
    pp_rot[2]=-pp_daugters[0]*sin(theta_in)+pp_daugters[2]*cos(theta_in);
    pp_rot[5]=-pp_daugters[3]*sin(theta_in)+pp_daugters[5]*cos(theta_in);

    double dif = pin[0]-pp_rot[0]-pp_rot[3]+pin[1]-pp_rot[1]-pp_rot[4]+pin[2]-pp_rot[2]-pp_rot[5];
    if(dif>0.001) cout<<"FUCK: "<<pin[0]<<","<<pin[1]<<","<<pin[2]<<","<<pp_rot[0]<<","<<pp_rot[1]<<","<<pp_rot[2]<<","<<pp_rot[3]<<","<<pp_rot[4]<<","<<pp_rot[5]<<endl;

    return pp_rot;

}
