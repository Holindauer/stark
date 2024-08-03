use crate::modules::field::{*};
use crate::modules::univariate_poly::{*};
use crate::modules::multivariate_poly::{*};
use num_bigint::RandBigInt;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use num_traits::ToPrimitive;

// Rescue Prime Hash Function
pub struct RescuePrime {
    pub p: BigInt,   // prime 
    pub m: usize,    // state width
    pub rate: usize, // rate for spone func
    pub capacity: usize,
    pub N: usize, 
    pub alpha: BigInt, // smallest invertible power
    pub alpha_inv: BigInt,
    pub MDS: Vec<Vec<FieldElement>>,
    pub MDS_inv: Vec<Vec<FieldElement>>,    
    pub round_constants: Vec<FieldElement>,
}


impl RescuePrime {

    // constructor
    pub fn new()-> RescuePrime {
        // setup parameters
        let p: BigInt = BigInt::from(407) * (BigInt::from(2).pow(119)) + 1;
        let m: usize = 2;
        let rate: usize = 1;
        let capacity: usize = 1;
        let N: usize = 27;
        let alpha: BigInt = BigInt::from(3);
        let alpha_inv: BigInt = BigInt::from(180331931428153586757283157844700080811 as u128);

        // MDS matrix
        let MDS: Vec<Vec<FieldElement>> = vec![
            vec![FieldElement::new(BigInt::from(270497897142230380135924736767050121214 as u128)), FieldElement::new(BigInt::from(4))],
            vec![FieldElement::new(BigInt::from(270497897142230380135924736767050121205 as u128)), FieldElement::new(BigInt::from(13))],
        ];

        // MDS inverse matrix
        let MDS_inv: Vec<Vec<FieldElement>> = vec![
            vec![FieldElement::new(BigInt::from(210387253332845851216830350818816760948 as u128)), FieldElement::new(BigInt::from(60110643809384528919094385948233360270 as u128))],
            vec![FieldElement::new(BigInt::from(90165965714076793378641578922350040407 as u128)), FieldElement::new(BigInt::from(180331931428153586757283157844700080811 as u128))],
        ];

        // pre field element round constants
        let round_constant_values: Vec<u128> = vec![
            174420698556543096520990950387834928928,
            109797589356993153279775383318666383471,
            228209559001143551442223248324541026000,
            268065703411175077628483247596226793933,
            250145786294793103303712876509736552288,
            154077925986488943960463842753819802236,
            204351119916823989032262966063401835731,
            57645879694647124999765652767459586992,
            102595110702094480597072290517349480965,
            8547439040206095323896524760274454544,
            50572190394727023982626065566525285390,
            87212354645973284136664042673979287772,
            64194686442324278631544434661927384193,
            23568247650578792137833165499572533289,
            264007385962234849237916966106429729444,
            227358300354534643391164539784212796168,
            179708233992972292788270914486717436725,
            102544935062767739638603684272741145148,
            65916940568893052493361867756647855734,
            144640159807528060664543800548526463356,
            58854991566939066418297427463486407598,
            144030533171309201969715569323510469388,
            264508722432906572066373216583268225708,
            22822825100935314666408731317941213728,
            33847779135505989201180138242500409760,
            146019284593100673590036640208621384175,
            51518045467620803302456472369449375741,
            73980612169525564135758195254813968438,
            31385101081646507577789564023348734881,
            270440021758749482599657914695597186347,
            185230877992845332344172234234093900282,
            210581925261995303483700331833844461519,
            233206235520000865382510460029939548462,
            178264060478215643105832556466392228683,
            69838834175855952450551936238929375468,
            75130152423898813192534713014890860884,
            59548275327570508231574439445023390415,
            43940979610564284967906719248029560342,
            95698099945510403318638730212513975543,
            77477281413246683919638580088082585351,
            206782304337497407273753387483545866988,
            141354674678885463410629926929791411677,
            19199940390616847185791261689448703536,
            177613618019817222931832611307175416361,
            267907751104005095811361156810067173120,
            33296937002574626161968730356414562829,
            63869971087730263431297345514089710163,
            200481282361858638356211874793723910968,
            69328322389827264175963301685224506573,
            239701591437699235962505536113880102063,
            17960711445525398132996203513667829940,
            219475635972825920849300179026969104558,
            230038611061931950901316413728344422823,
            149446814906994196814403811767389273580,
            25535582028106779796087284957910475912,
            93289417880348777872263904150910422367,
            4779480286211196984451238384230810357,
            208762241641328369347598009494500117007,
            34228805619823025763071411313049761059,
            158261639460060679368122984607245246072,
            65048656051037025727800046057154042857,
            134082885477766198947293095565706395050,
            23967684755547703714152865513907888630,
            8509910504689758897218307536423349149,
            232305018091414643115319608123377855094,
            170072389454430682177687789261779760420,
            62135161769871915508973643543011377095,
            15206455074148527786017895403501783555,
            201789266626211748844060539344508876901,
            179184798347291033565902633932801007181,
            9615415305648972863990712807943643216,
            95833504353120759807903032286346974132,
            181975981662825791627439958531194157276,
            267590267548392311337348990085222348350,
            49899900194200760923895805362651210299,
            89154519171560176870922732825690870368,
            265649728290587561988835145059696796797,
            140583850659111280842212115981043548773,
            266613908274746297875734026718148328473,
            236645120614796645424209995934912005038,
            265994065390091692951198742962775551587,
            59082836245981276360468435361137847418,
            26520064393601763202002257967586372271,
            108781692876845940775123575518154991932,
            138658034947980464912436420092172339656,
            45127926643030464660360100330441456786,
            210648707238405606524318597107528368459,
            42375307814689058540930810881506327698,
            237653383836912953043082350232373669114,
            236638771475482562810484106048928039069,
            168366677297979943348866069441526047857,
            195301262267610361172900534545341678525,
            2123819604855435621395010720102555908,
            96986567016099155020743003059932893278,
            248057324456138589201107100302767574618,
            198550227406618432920989444844179399959,
            177812676254201468976352471992022853250,
            211374136170376198628213577084029234846,
            105785712445518775732830634260671010540,
            122179368175793934687780753063673096166,
            126848216361173160497844444214866193172,
            22264167580742653700039698161547403113,
            234275908658634858929918842923795514466,
            189409811294589697028796856023159619258,
            75017033107075630953974011872571911999,
            144945344860351075586575129489570116296,
            261991152616933455169437121254310265934,
            18450316039330448878816627264054416127
        ];

        // convert round constant values to field elements
        let mut round_constants: Vec<FieldElement> = Vec::new();
        for constant in round_constant_values {
            let constant_field_element = FieldElement::new(BigInt::from(constant));
            round_constants.push(constant_field_element);
        }

        RescuePrime {p, m, rate, capacity, N, alpha, alpha_inv, MDS, MDS_inv, round_constants,}
    }

    // hash
    pub fn hash(&self, input_element: FieldElement) -> FieldElement {

        // absorb
        let mut state: Vec<FieldElement> = Vec::new();
        state.push(input_element);
        for _ in 0..(self.m-1) { state.push(FieldElement::zero()); }

        // permutation
        for r in 0..self.N {
            
            // forward half-round
            // S-box
            for i in 0..self.m {
                state[i] = state[i].pow(self.alpha.clone().to_u128().unwrap());
            }   

            // matrix
            let mut temp: Vec<FieldElement> = vec![];
            for _ in 0..self.m { temp.push(FieldElement::zero()); }
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] = temp[i].clone() + self.MDS[i][j].clone() * state[j].clone();
                }
            }

            // constants
            state = vec![];
            for i in 0..self.m {
                state.push( temp[i].clone() + self.round_constants[2*r*self.m+i].clone() );
            }

            // backward half-round
            // S-box
            for i in 0..self.m {
                state[i] = state[i].pow(self.alpha_inv.clone().to_u128().unwrap());
            }

            // matrix 
            temp = vec![];
            for _ in 0..self.m { temp.push(FieldElement::zero()); }
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] = temp[i].clone() + self.MDS[i][j].clone() * state[j].clone(); 
                }
            }

            // constants
            state = vec![];
            for i in 0..self.m {
                state.push( temp[i].clone() + self.round_constants[2*r*self.m+self.m+i].clone() );
            }
        }

        // squeeze
        state[0].clone()
    }       

    // trace
    pub fn trace(&self, input_element: FieldElement) -> Vec<Vec<FieldElement>> {

        // init trace
        let mut trace: Vec<Vec<FieldElement>> = vec![];

        // absorb
        let mut state: Vec<FieldElement> = Vec::new();
        state.push(input_element);
        for _ in 0..(self.m-1) { state.push(FieldElement::zero()); }

        // explicit copy to record state into trace
        trace.push(state.clone());

        // permutation
        for r in 0..self.N {

            // forward half-round
            // S-box
            for i in 0..self.m {
                state[i] = state[i].pow(self.alpha.clone().to_u128().unwrap());
            }

            // matrix
            let mut temp: Vec<FieldElement> = vec![];
            for _ in 0..self.m { temp.push(FieldElement::zero()); }
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] = temp[i].clone() + self.MDS[i][j].clone() * state[j].clone();
                }
            }

            // constants
            state = vec![];
            for i in 0..self.m {
                state.push( temp[i].clone() + self.round_constants[2*r*self.m+i].clone() );
            }

            // backward half-round
            // S-box
            for i in 0..self.m {
                state[i] = state[i].pow(self.alpha_inv.clone().to_u128().unwrap());
            }

            // matrix 
            temp = vec![];
            for _ in 0..self.m { temp.push(FieldElement::zero()); }
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] = temp[i].clone() + self.MDS[i][j].clone() * state[j].clone(); 
                }
            }

            // constants
            state = vec![];
            for i in 0..self.m {
                state.push( temp[i].clone() + self.round_constants[2*r*self.m+self.m+i].clone() );
            }

            // record state into trace
            trace.push(state.clone());
        }

        trace
    }

    // boundary constraints
    pub fn boundary_constraints(&self, output_element: FieldElement) -> Vec<(usize, usize, FieldElement)> {
        let mut constraints: Vec<(usize, usize, FieldElement)> = vec![];
        
        // at stark, capacity is zero
        constraints.push( (0, 1, FieldElement::zero()) );
        
        // at end, rate part is the fiven output element
        constraints.push( (self.N, 0, output_element) );

        constraints
    }

    // round constant polynomials
    pub fn round_constants_polynomials(&self, omicron: FieldElement) -> (Vec<MPolynomial>, Vec<MPolynomial>) {

        // ! Issue here

        // first step constants
        let mut first_step_constants: Vec<MPolynomial> = vec![];
        for i in 0..self.m {    

            // domain
            let mut domain: Vec<FieldElement> = vec![];
            for r in 0..self.N { domain.push(omicron.pow(r as u128)); }

            // values
            let mut values: Vec<FieldElement> = vec![];
            for r in 0..self.N { values.push(self.round_constants[2*r*self.m+i].clone()); }

            // interpolate univariate poly over domain and values
            let univariate: Polynomial = Polynomial::lagrange(domain, values);

            // lift uni to multivariate
            let multivariate: MPolynomial = MPolynomial::lift(&univariate, 0);


            // push to first step constants
            first_step_constants.push(multivariate);
        }

        // second step constants
        let mut second_step_constants: Vec<MPolynomial> = vec![];
        for i in 0..self.m {

            // domain
            let mut domain: Vec<FieldElement> = vec![];
            for r in 0..self.N { domain.push( omicron.pow(r as u128)); }

            // values
            let mut values: Vec<FieldElement> = vec![];
            for r in 0..self.N { values.push(self.round_constants[2*r*self.m+self.m+i].clone()); }

            // interpolate uni poly over domain
            let univariate: Polynomial = Polynomial::lagrange(domain, values);

            // lift uni to multivariate
            let multivariate: MPolynomial = MPolynomial::lift(&univariate, 0);

            // push to second step constants
            second_step_constants.push(multivariate);
        }

        (first_step_constants, second_step_constants)
    }


    // transition constraints
    pub fn transition_constraints(&self, omicron: FieldElement) -> Vec<MPolynomial> {
        // get polynomails that interpolate through the round constants
        let (first_step_constants, second_step_constants) = self.round_constants_polynomials(omicron);

        // arithmetize one round of Rescue-Prime
        let variables: Vec<MPolynomial> = MPolynomial::variables(1 + 2*self.m);
        let previous_state = &variables[1..(1+self.m)];
        let next_state = &variables[(1+self.m)..(1+2*self.m)];
        let mut air: Vec<MPolynomial> = vec![];

        for i in 0..self.m {

            // compute left hand side symbolically
            let mut lhs = MPolynomial::constant(0);
            for k in 0..self.m {
                lhs = lhs + MPolynomial::constant(self.MDS[i][k].value.to_u128().unwrap()) * (previous_state[k].clone().pow(self.alpha.to_u128().unwrap()))
            }
            lhs = lhs + first_step_constants[i].clone();

            // compute right hand side symbolically
            let mut rhs = MPolynomial::constant(0);
            for k in 0..self.m {
                rhs = rhs + MPolynomial::constant(self.MDS_inv[i][k].value.to_u128().unwrap()) * (next_state[k].clone() - second_step_constants[k].clone());
            }
            rhs = rhs.pow(self.alpha.to_u128().unwrap());

            // equate left and right hand sides
            air.push(lhs - rhs); // should equal zero
        }

        air
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::field::{*};
    use crate::modules::univariate_poly::{*};
    use crate::modules::multivariate_poly::{*};
    use num_bigint::RandBigInt;
    use num_bigint::BigInt;
    use num_traits::{Zero, One};
    use rand::RngCore;
    use rand::Rng;
    use rand::random;

    #[test]
    fn test_rescue_prime() {
        let rp = RescuePrime::new();

        // test vectors
        assert_eq!(
            rp.hash(FieldElement::one()),
            FieldElement::new(BigInt::from(244180265933090377212304188905974087294 as u128))
        );
        assert_eq!(
            rp.hash(FieldElement::new(BigInt::from(57322816861100832358702415967512842988 as u128))),
            FieldElement::new(BigInt::from(89633745865384635541695204788332415101 as u128))
        );

        // test trace boundaries
        let a = FieldElement::new(BigInt::from(57322816861100832358702415967512842988 as u128));
        let b = FieldElement::new(BigInt::from(89633745865384635541695204788332415101 as u128));
        let trace = rp.trace(a.clone());

        assert_eq!(trace[0][0], a);
        assert_eq!(trace[trace.len()-1][0], b);
    }   

    #[test]
    fn test_trace() {

        // init rescue prime
        let rp = RescuePrime::new();

        // hash test
        let input_element = FieldElement::new(BigInt::from(57322816861100832358702415967512842988 as u128));
        let b = FieldElement::new(BigInt::from(89633745865384635541695204788332415101 as u128));
        let output_element = rp.hash(input_element.clone());
        assert_eq!(b, output_element);
        
        // get trace
        let mut trace = rp.trace(input_element);

        // test boundary constraints
        let constraints = rp.boundary_constraints(output_element.clone());
        for condition in constraints {
            let (cycle, element, value) = condition;
            if trace[cycle][element] != value {
                println!("Rescue prime boundary condition error: trace element: {} at cycle {} has value {}, but should have value {}", element, cycle, trace[cycle][element], value);
                assert!(false);
            }
        }

        // test transition constraints
        let omicron = FieldElement::primitive_nth_root(1 << 119);
        let transition_constraints: Vec<MPolynomial> = rp.transition_constraints(omicron.clone());
        let (first_step_constants, second_step_constants) = rp.round_constants_polynomials(omicron.clone());
        for o in 0..trace.len()-1 {
            println!("cycle: {}", o);
            for air_poly in rp.transition_constraints(omicron.clone()) {
                
                // get prev and next state
                let previous_state: Vec<FieldElement> = vec![trace[o][0].clone(), trace[o][1].clone()];
                let next_state: Vec<FieldElement> = vec![trace[o+1][0].clone(), trace[o+1][1].clone()];
                
                // get point
                let mut point = vec![omicron.clone().pow(o as u128)];
                point.extend(previous_state.clone());
                point.extend(next_state.clone());

                println!("air poly {:?}", air_poly);

                if air_poly.eval(&point) != FieldElement::zero() {
                    println!("Rescue prime transition condition error: air polynomial does not evaluate to zero at point: {:?}", point);
                    assert!(false);
                }
            }
        }   

        // insert errors into trace to make sure errors get notices
        for k in 0..10 { 
            println!("trial: {}...", k);

            // sample error location and value randomly
            let mut register_index: usize = random();
            register_index = register_index % rp.m;
            let mut cycle_index: usize = random();
            cycle_index = cycle_index % (rp.N+1);

            // random value
            let mut rng = rand::thread_rng();
            let mut random_bytes: Vec<u8> = vec![0; 17];
            rng.fill_bytes(&mut random_bytes);
            let mut value: FieldElement = FieldElement::sample(random_bytes);

            // skip if value is zero
            if value == FieldElement::zero() {
                continue;
            }

            // reproduce deterministic error 
            if k == 0 {
                register_index = 1;
                cycle_index = 22;
                value = FieldElement::new(BigInt::from(17274817952119230544216945715808633996 as u128));
            }
            
            // perturb
            trace[cycle_index][register_index] = trace[cycle_index][register_index].clone() + value.clone();

            // flag for error noticed
            let mut error_got_noticed: bool = false;

            // test boudnary constraints
            println!("testing boundary constraints...");
            for condition in rp.boundary_constraints(output_element.clone()) {
                if error_got_noticed { break; }

                // get constraint
                let (cycle, element, value) = condition;

                // if trace doesn't match, flag error
                if trace[cycle][element] != value {
                    error_got_noticed = true;
                }
            }

            // test transition constraints
            println!("testing transition constraints...");
            for o in 0..trace.len()-1 {
                println!("cycle: {}", o);
                if error_got_noticed { break; }
                
                for air_poly in rp.transition_constraints(omicron.clone()){

                    // get next and prev state
                    let previous_state = vec![trace[o][0].clone(), trace[o][1].clone()];
                    let next_state = vec![trace[o+1][0].clone(), trace[o+1][1].clone()];  

                    // get point
                    let mut point = vec![omicron.clone().pow(o as u128)];
                    point.extend(previous_state.clone());
                    point.extend(next_state.clone());  

                    // if air polynomial doesn't evaluate to zero, flag error
                    if air_poly.eval(&point) != FieldElement::zero() {
                        error_got_noticed = true;
                    }
                }
            }

            // if error got noticed, panic
            if !error_got_noticed {
                println!("error not noticed");
                println!("cycle: {}", cycle_index.clone());
                println!("register: {}", register_index.clone());
                println!("value: {}", value.clone());
                assert!(false);
            }

            // return trace to original state
            trace[cycle_index][register_index] = trace[cycle_index][register_index].clone() - value.clone();
        }

        println!("rescue prime test passed");

    }
}