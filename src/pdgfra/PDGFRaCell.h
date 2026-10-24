// PDGFRaCell.h : thin headless wrapper around the PDGFa-project CCellPart
// (mouse ICC-derived cell with ANO1/CaL/VDDR off, stochastic SK + IP3R),
// used as a PDGFRa+ cell. Pimpl: this header does NOT include CellPart.h,
// so it can be mixed with headers that define twoelement/fourelement/sixelement.
#pragma once
#include <memory>

class CCellPart;

class PDGFRaCell
{
public:
	explicit PDGFRaCell(unsigned seed = 1);
	PDGFRaCell(const PDGFRaCell& o);             // deep copy (independent cell, same RNG state)
	PDGFRaCell& operator=(const PDGFRaCell& o);  // deep copy
	PDGFRaCell(PDGFRaCell&&) noexcept;
	PDGFRaCell& operator=(PDGFRaCell&&) noexcept;
	~PDGFRaCell();

	void Reseed(unsigned seed);  // reseed this cell's private RNG
	unsigned Seed() const;

	// InitialCellVariables() + Nae/Ke/Cae = 140/5.4/1.8 mM + stochastic channel
	// initialisation (as in the MFC driver) + timeX = 0.
	void Init();

	// Advance exactly dt (ms) in current-clamp mode. Iext_pA is outward-positive
	// external current (pA) added to the membrane equation dV/dt = -(Itot+Iext)/Cap.
	// If a step produces a negative concentration it is retried with dt*0.1
	// sub-steps (as the GUI driver does) until the full dt has been covered.
	// Returns false if the step could not be completed (sub-step < 1e-9 ms).
	bool Step(double dt, double Iext_pA);

	double Vm() const;     // mV
	double Cai() const;    // mM (bulk cytosol free Ca)
	double ISK() const;    // pA (outward positive)
	double IP3() const;    // mM
	double Cap() const;    // pF
	double Time() const;   // ms (cell's internal timeX)
	double CaPMU() const;  // mM (pacemaker-unit Ca)
	double Itot() const;   // pA, total ionic current of the last step
	long long Retries() const;  // cumulative number of dt*0.1 retries

	// neural hook: PIP2->IP3 production rate multiplied by (1+s); default 0
	void SetAgonist(double s);
	double Agonist() const;
	// SK conductance scale factor (SF_SK, default 1)
	void SetSKScale(double s);
	// [modified PDGFa] SK reversal = EK, SK voltage-independent, and scale factors for background inward currents
	void SetModified(bool skUseEK, bool skVoltIndep, double sfNaLCN, double sfSOC);
	// [modified PDGFa] purinergic drive: target agonist-evoked IP3 (mM) added to the IP3 seen by IP3R (tau 300 ms)
	void SetIP3Agonist(double mM);
	// [modified PDGFa] SK Ca2+ source: 0 original SOC nanodomain, 1 submembrane PMU Ca2+ with Kd (mM) and Hill coefficient
	void SetSKCaSource(int source, double kd_mM, double hill);

	// escape hatch for advanced use (include CellPart.h to use it)
	CCellPart& Raw();
	const CCellPart& Raw() const;

private:
	std::unique_ptr<CCellPart> p;
	long long nRetries = 0;
};
