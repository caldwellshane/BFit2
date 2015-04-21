{
	TFile *tfile = new TFile("BFit.root");
	
	Int_t test	= 0;
	Int_t data	= 1;
	Int_t MC	= 0;
	
	if (test) {
		cTest->Close(); cTest->Draw();
	}
	if (data) {
		c_BFit->Close(); c_BFit->Draw();
	}
	if (MC) {
//		c_BFit->Close(); c_BFit->Draw();
//		c_feeding->Close(); c_feeding->Draw();
//		c_decays_cyctime->Close(); c_decays_cyctime->Draw();
		c_Ti->Close(); c_Ti->Draw();
		c_Vi->Close(); c_Vi->Draw();
		c_Wi->Close(); c_Wi->Draw();
		c_Zi->Close(); c_Zi->Draw();
		c_Xi->Close(); c_Xi->Draw();
		c_Yi->Close(); c_Yi->Draw();
//		c_Ui->Close(); c_Ui->Draw();
//		c_T1V1W1Z1X2->Close(); c_T1V1W1Z1X2->Draw();
//		c_T2V2W2Z2X3Y3->Close(); c_T2V2W2Z2X3Y3->Draw();
//		c_T3V3W3Z3->Close(); c_T3V3W3Z3->Draw();
	}
}
