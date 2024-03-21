# Application submission info

---
> Submitted on 3-Mar-2024; approved 5-Mar-2024 (BIO240064)
---

Submitted through ACCESS: https://access-ci.org/

Submitted an "Explore ACCESS", different types listed here: https://allocations.access-ci.org/project-types

After logged in, submission process started from this page: https://allocations.access-ci.org/opportunities

Then selecting to "Request New Project", and choosing "Request an Explore ACCESS project".

## Required information I added to the form

**Title:** STAMPS 2024 at Marine Biologial Laboratory in Woods Hole, MA, USA

**Public overview:**  

> General overview
> The STAMPS course (Strategies and Techniques for Analyzing Microbial Population Structures) has been a yearly event at the Marine Biological Laboratory (MBL) in Woods Hole, MA, USA for over a decade, helping nearly 1,000 learners establish a foundation in bioinformatics over the years. Since 2019 we have been fortunate enough to have our computational infrastructure provided through XSEDE, then ACCESS, and Jetstream (1 and 2), and we hope to continue that this year. The course takes place this year 17-July to 28-July and will involve a lot of command-line, R, amplicon analysis, metagenomics, statistical training, and other concepts and computational activities (last year's schedule can be seen here: https://github.com/mblstamps/stamps2023/wiki#schedule).
> 
> How we plan to use ACCESS resources
> We will be running an about 11-day workshop (19-July to 29-July). I have experience in the past doing this through Indiana Jetstream2, managing through Exosphere, so that would be great if possible. We are requesting the ability to run 60 (accounting for faculty and participants) m3.large instances concurrently for about 15 days straight. This is longer than we plan on them being active, but would like to leave room for building the image and testing/setting things up ahead of time, and any unexpected situations. Using the Jetstream2 estimator (https://docs.jetstream-cloud.org/alloc/estimator/), with 60 m3.large instances for 15 days, this comes out to a request for 345,600 SUs/ACCESS Credits. We think this will work under the 400,000 max for "Explore ACCESS" allocations. Since we will be using these all at once over a short time, we do request the full amount to be released up front if possible. Please let me know if I should be estimating a different way or if any other information is needed.
> 
> Thanks for your consideration and any help!

**Keywords:** microbial ecology, bioinformatics

**Grant type:** other

**Opportunity questions:** didn’t check any (though these options seem to change over time, so see if any are relevant)

**Fields of Science:** “Other Biological Sciences”

**Documents:** Needed to attach my CV here.

**Available resources:** checked box for ACCESS Credits

---

## After submission

### Getting the remaining 200,000 credits
As of 2024, maybe 2023, the “Explore ACCESS” project comes with a total of 400,000 credits, with 200,000 released at first. After getting approved, I went to the https://allocations.access-ci.org/requests page, selected the STAMPS 2024 project, clicked on “Credits + Resources”, then “REQUEST MORE CREDITS”, then “REQUEST MORE CREDITS” again, then “REQUEST A SUPPLEMENT”. Then filled out these:
- Reason: 
  -	This is for a course that will be using all the credits over a short period of time, so we would like the total 400,000 available from the start. Thanks!
- Available Resources:
  -	Checked the box for ACCESS Credits
- Document Type:
  -	A Progress Report is required. This was just a PDF state the same as the “Reason” above.

Then submitted the form.

### Transferring credits from ACCESS to JetStream2
Once approved, and logged in, needed this page (https://allocations.access-ci.org/requests) in order to transfer ACCESS credits to a specific resource. For the appropriate allocation/Project, selected "Credits + Resources", then the text box that initially says "Add a resource to your exchange...", then selected "Indiana Jetstream2 CPU", then entered 196,000 (since I did this before the other 200,000 were available). Then for "Indiana Jetstream2 Storage" added the remaining 4,000, giving 4 TB of shared storage to use.
After the other 200,000 are made available, repeat the process without adding anymore to Storage.

### Requesting quota-limit increases so we can run up to 60 instances concurrently
The allocation comes with a limit on the number of concurrent instances that can be run. I submitted a request when logged into Jetstream2 here: https://jetstream2.exosphere.app/exosphere/getsupport

This is the text I submitted (after selecting the radio dial for “An Allocation”):

> Hi there :)
> 
> We plan to use this allocation (BIO240064) with 60 concurrent m3.large instances for a bioinformatics course we are running. The starting limit is set to 10.
> 
> Could you please help with increasing the allotted quotas so that we will be able to run up to 60 m3.large instances concurrently on this allocation, including cores, ram, volume, ports, available IP addresses, and whatever other magic you folks take care of that I'm naive to?
>
> Thank you for any help!  
> -Mike
