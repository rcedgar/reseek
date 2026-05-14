#include "myutils.h"

static const string default_varstr =
"aa20=5.29E-01;angle32=8.66E-02;dali=4.35E-04;gap2=8.94E-01;lddt=4.66E-02;mendist32=9.01E-03;nendist32=1.16E-01;nensec32=1.58E-02;pack32=8.87E-02;pendist32=1.37E-02;pm2=1.40E-02;pmdd32=8.04E-03;ppack32=1.77E-03;rendist32=2.25E-02;revw=6.72E-01;sec32=5.97E-02;selfw=9.62E-01;turnd32=3.49E-02;";

void parse_varstr(
	const string &arg_VarStr,
	vector<string> &Names,
	vector<float> &Values)
	{
	Names.clear();
	Values.clear();

	const string &VarStr = (arg_VarStr == "" ? default_varstr : arg_VarStr);

	vector<string> Fields;
	Split(VarStr, Fields, ';');

	const uint n = SIZE(Fields);
	for (uint i = 0; i < n; ++i)
		{
		const string &NameEqValue = Fields[i];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("SubsetBench::parse_varstr(%s) not name=value '%s'",
				VarStr.c_str(), Fields[i].c_str());
		const string &Name = Fields2[0];
		const string &ValueStr = Fields2[1];
		float Weight = StrToFloatf(ValueStr);
		Names.push_back(Name);
		Values.push_back(Weight);
		}
	}
