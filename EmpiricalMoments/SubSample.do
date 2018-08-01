**************************ENAHO*********************************************
{
clear all
set more off
cd "C:\Users\razuero\Dropbox\OptmalTaxationShared\Data\DataAnalysis\All\ENAHO\ENAHOARM"
use PER_2007a_BID.dta

*-----------------------*
*Organizing the dataset *
*-----------------------*

*Keep Lima
keep if region_c==15

*Trim top 1%. 
*Si es ocupado
sum ylmpri_ci if emp_ci!=0, d
drop if ylmpri_ci>r(p99) & ylmpri_ci!=.

/*Information about informality:
p599> trabajador independiente formal
p511a tipo de contrato

*Definición comunmente utilizada de informalidad:
-Patronos y cuenta propia cuya unidad productiva es informal
-Asalariados sin seguridad social financiada por empleador o sin contrato
-Trabajadores familiares no remunerados
*/


*Definición de informalidad. Si tienen contrato
g formal_empleado=1 if emp_ci==1 & categopri_ci==3 & (tipocontrato==1 | tipocontrato==2) 
replace formal_empleado=0 if emp_ci==1 & categopri_ci==3 & tipocontrato==3

*Número de personas que laboran en centro de trabajo
g numero_trabajadores=p512b
replace numero=. if numero== 9998 

export delimited using "C:\Users\razuero\Dropbox\OptmalTaxationShared\Data\DataAnalysis\All\ENAHO\ENAHOARM\OptimaltaxationSubSampleENAHO.csv", replace

}

**************************CENSO*********************************************
{
clear all
set more off
cd "C:\Users\razuero\Dropbox\OptmalTaxationShared\Data\DataAnalysis\All\Census\Modified"
use MergeCenso.dta

*--------------------------------------------------*
*Generating different ways of corporate income tax-*
*--------------------------------------------------*

*1. All corporate income tax that is reported directly
g CITAX=CAP5MONTO36 
order CITAX

*2. Replacing reports from the RER

*we assume that some people didn;t report paying corporate income tax
*when they payed RER. We will take it as corporate income TAX. 
replace CITAX=RER3 if RER3!=. & RER3!=0 & RER3>CITAX

*3. RUS
*We also consider corporate income tax the payments that are given to the RUS
count if CAP5MONTO36!=0 & CAP5MONTO36!=. & RUS3!=0 & RUS3!=.
replace CITAX=RUS3 if RUS3!=0 & RUS3!=.

*4. CITAX in USD
g CITAXUSD=CITAX*0.315



*================*
*Sample selection*
*================*

/* - - - - - - - - - - - */

keep if PROVINCIA=="LIMA"


*Trimming of PROFITS*
sum CITAX, d
drop if CITAX> r(p99)

*Eliminating firms that report negative taxes
drop if CITAX<0

*Excluding the ones with no revenue in 2007. 
g TOTINCOMEFORMAL=CAP5MONTO1/12+CAP5MONTO5/12
drop if TOTINCOMEFORMAL<=0

*Exclude firms that report not having workers
drop if NTRABAJADORES8==0

*Exclude firms without production
drop if CAP5MONTO9<=0
*Only 5 firms dropped

export delimited using "C:\Users\razuero\Dropbox\OptmalTaxationShared\Data\DataAnalysis\All\Census\Modified\OptimaltaxationSubSampleCenso.csv", replace

}
