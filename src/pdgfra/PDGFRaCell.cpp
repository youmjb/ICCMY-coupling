// PDGFRaCell.cpp : see PDGFRaCell.h
#include "PDGFRaCell.h"
#include "CellPart.h"
#include <cmath>

PDGFRaCell::PDGFRaCell(unsigned seed) : p(new CCellPart(seed)) {}
PDGFRaCell::PDGFRaCell(const PDGFRaCell& o) : p(new CCellPart(*o.p)), nRetries(o.nRetries) {}
PDGFRaCell& PDGFRaCell::operator=(const PDGFRaCell& o)
{
	if (this != &o) { *p = *o.p; nRetries = o.nRetries; }
	return *this;
}
PDGFRaCell::PDGFRaCell(PDGFRaCell&&) noexcept = default;
PDGFRaCell& PDGFRaCell::operator=(PDGFRaCell&&) noexcept = default;
PDGFRaCell::~PDGFRaCell() = default;

void PDGFRaCell::Reseed(unsigned seed) { p->Reseed(seed); }
unsigned PDGFRaCell::Seed() const { return p->seed; }

void PDGFRaCell::Init()
{
	CCellPart& c = *p;
	c.InitialCellVariables();          // also calls InitializeStochasticChannels()
	c.Nae = 140.0;                      // set by the UI in the MFC program
	c.Ke = 5.4;
	c.Cae = 1.8;
	c.InitializeStochasticChannels();   // MitoView calls this again before the run loop
	c.timeX = 0.0;
	c.timeX0 = 0.0;
	nRetries = 0;
}

bool PDGFRaCell::Step(double dt, double Iext_pA)
{
	CCellPart& c = *p;
	const double t0 = c.timeX;
	const double tEnd = t0 + dt;
	double remaining = dt;
	const double minDt = 1e-9;
	while (remaining > 1e-12 * dt)
	{
		c.m_dDt = remaining;
		while (!c.RunCellModel(true, Iext_pA))
		{
			c.m_dDt *= 0.1;
			++nRetries;
			if (c.m_dDt < minDt)
			{
				c.timeX = tEnd;   // keep the time axis consistent; state held
				return false;
			}
		}
		c.timeX += c.m_dDt;
		remaining -= c.m_dDt;
	}
	c.timeX = tEnd;   // remove round-off of the sub-step sum
	return true;
}

double PDGFRaCell::Vm() const { return p->Vm; }
double PDGFRaCell::Cai() const { return p->Cai; }
double PDGFRaCell::ISK() const { return p->ISK; }
double PDGFRaCell::IP3() const { return p->IP3; }
double PDGFRaCell::Cap() const { return ::Cap; }
double PDGFRaCell::Time() const { return p->timeX; }
double PDGFRaCell::CaPMU() const { return p->CaPMU; }
double PDGFRaCell::Itot() const { return p->Itot; }
long long PDGFRaCell::Retries() const { return nRetries; }

void PDGFRaCell::SetAgonist(double s) { p->agonistScale = s; }
double PDGFRaCell::Agonist() const { return p->agonistScale; }
void PDGFRaCell::SetSKScale(double s) { p->SF_SK = s; }
void PDGFRaCell::SetModified(bool ek, bool vi, double sfNaLCN, double sfSOC) { p->skUseEK = ek; p->skVoltIndep = vi; p->gapCarriedByK = ek; p->skDeterministic = ek; p->SF_NaLCN = sfNaLCN; p->SF_SOC = sfSOC; }

CCellPart& PDGFRaCell::Raw() { return *p; }
const CCellPart& PDGFRaCell::Raw() const { return *p; }
void PDGFRaCell::SetIP3Agonist(double mM) { p->ip3AgTarget = mM; }
void PDGFRaCell::SetSKCaSource(int source, double kd, double hill) { p->skCaSource = source; p->skKd = kd; p->skHill = hill; }
