{
TChain chain("tree") ;
chain.Add("scattering_events_17F_284.root");
chain.GetListOfFiles()->Print();
chain.MakeClass("AnaClass2");
}
