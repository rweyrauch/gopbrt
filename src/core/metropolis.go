package core

import()

type (
	
   	MetropolisRenderer struct {
        camera Camera
        bidirectional bool
        nDirectPixelSamples, nPixelSamples, maxDepth int
        largeStepsPerPixel, nBootstrap, maxConsecutiveRejects int
        directLighting *DirectLightingIntegrator
        nTasksFinished int       
   }

)

func NewMetropolisRenderer(perPixelSamples, nBootstrap, nDirectPixelSamples int, largeStepProbability float64, doDirectSeparately bool,
        maxrejects, maxdepth int, camera Camera, doBidirectional bool) *MetropolisRenderer {
	renderer := new(MetropolisRenderer)
	renderer.camera = camera
	renderer.bidirectional = doBidirectional
	renderer.nDirectPixelSamples = nDirectPixelSamples
	renderer.nPixelSamples = perPixelSamples
	renderer.maxDepth = maxdepth
	renderer.nBootstrap = nBootstrap
	renderer.maxConsecutiveRejects = maxrejects
		  
    renderer.largeStepsPerPixel = Maxi(1, int(RoundUpPow2(uint32(largeStepProbability * float64(renderer.nPixelSamples)))))
    if renderer.largeStepsPerPixel >= renderer.nPixelSamples {
    	 renderer.largeStepsPerPixel /= 2
    }	 
    //Assert(largeStepsPerPixel >= 1 && largeStepsPerPixel < nPixelSamples);
    if (renderer.nPixelSamples % renderer.largeStepsPerPixel) != 0 {
        origPixelSamples := renderer.nPixelSamples
        renderer.nPixelSamples += renderer.largeStepsPerPixel - (renderer.nPixelSamples % renderer.largeStepsPerPixel)
        Warning("Rounding up to %d Metropolis samples per pixel (from %d)",
                renderer.nPixelSamples, origPixelSamples)
    }

	if doDirectSeparately {	  
		//renderer.directLighting = NewDirectLightingIntegrator(SAMPLE_ALL_UNIFORM, renderer.maxDepth)
	}
	
	return renderer	
}
		
func (r *MetropolisRenderer) Render(scene *Scene) {
    //PBRT_MLT_STARTED_RENDERING();
    if len(scene.lights) > 0 {
        x0, x1, y0, y1 := r.camera.Film().GetPixelExtent()
        t0, t1 := r.camera.ShutterOpen(), r.camera.ShutterClose()
        lightDistribution := ComputeLightSamplingCDF(scene)

        if r.directLighting != nil {
            //PBRT_MLT_STARTED_DIRECTLIGHTING();
            // Compute direct lighting before Metropolis light transport
            if r.nDirectPixelSamples > 0 {        	
                sampler := NewLDSampler(x0, x1, y0, y1, r.nDirectPixelSamples, t0, t1)
                /*
                Sample *sample = new Sample(&sampler, directLighting, NULL, scene);
                vector<Task *> directTasks;
                int nDirectTasks = max(32 * NumSystemCores(),
                                 (r.camera.film.XResolution() * r.camera.film.YResolution()) / (16*16))
                nDirectTasks = RoundUpPow2(nDirectTasks)
                ProgressReporter directProgress(nDirectTasks, "Direct Lighting")
                for i := 0; i < nDirectTasks; i++ {
                    directTasks.push_back(new SamplerRendererTask(scene, this, camera, directProgress,
                                                                  &sampler, sample, false, i, nDirectTasks))
				}                    
                std::reverse(directTasks.begin(), directTasks.end())
                EnqueueTasks(directTasks)
                WaitForAllTasks()
                directProgress.Done()
				*/                
            }
            r.camera.Film().WriteImage(1.0)
            //PBRT_MLT_FINISHED_DIRECTLIGHTING();
        }
        // Take initial set of samples to compute $b$
        //PBRT_MLT_STARTED_BOOTSTRAPPING(nBootstrap);
        rng := CreateRNG(0)
        
        //var arena MemoryArena
        //vector<float> bootstrapI;
        //vector<PathVertex> cameraPath(maxDepth, PathVertex());
        //vector<PathVertex> lightPath(maxDepth, PathVertex());
        sumI := 0.0
        /*
        bootstrapI.reserve(nBootstrap);
        MLTSample sample(maxDepth);
        for i := 0; i < r.nBootstrap; i++ {
            // Generate random sample and path radiance for MLT bootstrapping
             x := Lerp(rng.RandomFloat(), x0, x1)
             y := Lerp(rng.RandomFloat(), y0, y1)
            r.LargeStep(rng, &sample, r.maxDepth, x, y, t0, t1, r.bidirectional)
            L := PathL(sample, scene, arena, r.camera, lightDistribution,
                               &cameraPath[0], &lightPath[0], rng)

            // Compute contribution for random sample for MLT bootstrapping
            float I = ::I(L);
            sumI += I
            bootstrapI.push_back(I);
            arena.FreeAll()
        }
        */
        b := sumI / float64(r.nBootstrap)
        //PBRT_MLT_FINISHED_BOOTSTRAPPING(b);
        Info("MLT computed b = %f", b);

        // Select initial sample from bootstrap samples
        contribOffset := rng.RandomFloat() * sumI
        rng.Seed(0)
        sumI = 0.0
        /*
        MLTSample initialSample(maxDepth);
        for i := 0; i < r.nBootstrap; i++ {
             x := Lerp(rng.RandomFloat(), x0, x1)
             y := Lerp(rng.RandomFloat(), y0, y1)
            r.LargeStep(rng, &initialSample, r.maxDepth, x, y, t0, t1, r.bidirectional)
            sumI += bootstrapI[i]
            if sumI > contribOffset {
                break
            }    
        }
		*/
        // Launch tasks to generate Metropolis samples
         nTasks := r.largeStepsPerPixel
        largeStepRate := r.nPixelSamples / r.largeStepsPerPixel
        Info("MLT running %d tasks, large step rate %d", nTasks, largeStepRate)
        /*
        ProgressReporter progress(nTasks * largeStepRate, "Metropolis");
        vector<Task *> tasks;
        Mutex *filmMutex = Mutex::Create();
        Assert(IsPowerOf2(nTasks));
        uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
        uint32_t pfreq = (x1-x0) * (y1-y0);
        for (uint32_t i = 0; i < nTasks; ++i) {
            float d[2];
            Sample02(i, scramble, d);
            tasks.push_back(new MLTTask(progress, pfreq, i,
                d[0], d[1], x0, x1, y0, y1, t0, t1, b, initialSample,
                scene, camera, this, filmMutex, lightDistribution));
        }
        EnqueueTasks(tasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < tasks.size(); ++i)
            delete tasks[i];
        progress.Done();
        Mutex::Destroy(filmMutex);
        delete lightDistribution;
		*/        
    }
    r.camera.Film().WriteImage(1.0)
    //PBRT_MLT_FINISHED_RENDERING();	
}

func (r *MetropolisRenderer) Li(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena, isect *Intersection, T *Spectrum) *Spectrum { return nil }
func (r *MetropolisRenderer) Transmittance(scene *Scene, ray *RayDifferential, sample *Sample, rng *RNG, arena *MemoryArena) *Spectrum { return nil }

func CreateMetropolisRenderer(params *ParamSet, camera Camera) *MetropolisRenderer { 
    largeStepProbability := params.FindFloatParam("largestepprobability", 0.25)
    perPixelSamples := params.FindIntParam("samplesperpixel", 100)
    nBootstrap := params.FindIntParam("bootstrapsamples", 100000)
    nDirectPixelSamples := params.FindIntParam("directsamples", 4)
    doDirectSeparately := params.FindBoolParam("dodirectseparately", true)
    mr := params.FindIntParam("maxconsecutiverejects", 512)
    md := params.FindIntParam("maxdepth", 7)
    doBidirectional := params.FindBoolParam("bidirectional", true)

    if options.QuickRender {
        perPixelSamples = Maxi(1, perPixelSamples / 4)
        nBootstrap = Maxi(1, nBootstrap / 4)
        nDirectPixelSamples = Maxi(1, nDirectPixelSamples / 4)
    }

    return NewMetropolisRenderer(perPixelSamples, nBootstrap,
        nDirectPixelSamples, largeStepProbability, doDirectSeparately,
        mr, md, camera, doBidirectional)
}
