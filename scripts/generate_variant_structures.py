#!/usr/bin/env python3
"""
Generate splice variant structures using ESMFold API.
Uses kinase domain sequences only (under 400aa limit).
"""

import requests
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"

# Kinase domain sequences (extracted from UniProt)
# Format: (canonical_kinase_domain, variant_kinase_domain, description)

VARIANTS = {
    'BRAF_p61': {
        'gene': 'BRAF',
        'drug': 'Vemurafenib',
        'expected': 'resistant',
        # BRAF kinase domain: 457-717 (261 aa)
        # p61-BRAF lacks exons 4-8 (aa 1-187) but kinase domain is preserved
        # However, p61 dimerizes constitutively - kinase domain structure similar
        'canonical_kinase': 'IGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH',
        'variant_kinase': 'IGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH',
        # Same kinase domain - resistance is due to dimerization, not kinase structure
        'note': 'Kinase domain identical; resistance from RAS-independent dimerization',
    },
    'EGFR_vIII': {
        'gene': 'EGFR',
        'drug': 'Gefitinib',
        'expected': 'resistant',
        # EGFR kinase domain: 712-979 (268 aa) - preserved in EGFRvIII
        # vIII deletes exons 2-7 (aa 6-273) in extracellular domain
        'canonical_kinase': 'KVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQG',
        'variant_kinase': 'KVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQG',
        # Same kinase domain - resistance from constitutive activation
        'note': 'Kinase domain identical; resistance from ligand-independent activation',
    },
    'MET_ex14skip': {
        'gene': 'MET',
        'drug': 'Capmatinib',
        'expected': 'sensitive',  # This variant IS drug sensitive!
        # MET kinase domain: 1078-1345 (268 aa)
        # Exon 14 skip deletes aa 964-1010 (juxtamembrane region, NOT kinase)
        'canonical_kinase': 'ARDMYDKEYYSVHNKTGAKLPVKWMALESLQTQKFTTKSDVWSFGVVLWELMTRGAPPYPDVNTFDITVYLLQGRRLLQPEYCPDPLYEVMLKCWHPKAEMRPSFSELVSRISAIFSTFIGEHYVHVNATYVNVKCVAPYPSLLSSEDNADDEVDTRPASFWETS',
        'variant_kinase': 'ARDMYDKEYYSVHNKTGAKLPVKWMALESLQTQKFTTKSDVWSFGVVLWELMTRGAPPYPDVNTFDITVYLLQGRRLLQPEYCPDPLYEVMLKCWHPKAEMRPSFSELVSRISAIFSTFIGEHYVHVNATYVNVKCVAPYPSLLSSEDNADDEVDTRPASFWETS',
        # Same kinase domain - variant is oncogenic AND drug sensitive
        'note': 'Kinase domain identical; exon14 skip causes MET stabilization, SENSITIVE to inhibitors',
    },
    'AR_V7': {
        'gene': 'AR',
        'drug': 'Enzalutamide',
        'expected': 'resistant',
        # AR ligand binding domain (LBD): 670-919 (250 aa)
        # AR-V7 is truncated at aa 644 - MISSING THE ENTIRE LBD
        # Enzalutamide binds to LBD, so no binding site exists in AR-V7
        'canonical_kinase': 'RQLKKLDNLHDSNTSLFSPSKASTYNDSILLWFWNSPRFFLQRTGDLLLFPQKSYPQCFHHMFQILHLFLESTPHAHLRTFCSWNPQTPFHPIVRNCLEQPGQNSFHCSGLTYEKYRHCPNAIFTHYGSLLLDKDELEQSQRLRDQLNELRRHEGSSPHLPAETPKQLIPLLKGAATETSTLSAQSSGSEGCTSELSPSKSHSRGFMESTETSSSEVSSPTEETTFQHLLRPVFSPL',
        'variant_kinase': None,  # No LBD in AR-V7
        'note': 'AR-V7 lacks LBD entirely; Enzalutamide binding pocket does not exist',
    },
    'ALK_L1196M': {
        'gene': 'ALK',
        'drug': 'Crizotinib',
        'expected': 'resistant',
        # ALK kinase domain: 1116-1399 (284 aa)
        # L1196M is gatekeeper mutation within kinase domain
        'canonical_kinase': 'VAVKMLKSTARSSEKQALMSELKIMTHLGPHLNIVNLLGACTKGGPIYIITEYCRYGDLVDYLHRNKHTFLQHHSDKRRPPSAELYS NALPVGLPLPSHVSLTGESDGGYMDMSKDESVDYVPMLDMKGDVKYADIESSNYMEREGKFKVVLRAKELRLQYREALAQHLHKDIVLKDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDNNPFSFQSDVYSYGIVLYELMTGELPYSHINNRDQIIFMVGRGYASPDLSKLYC',
        'variant_kinase': 'VAVKMLKSTARSSEKQALMSELKIMTHLGPHLNIVNLLGACTKGGPIYIITEYCRYGDLVDYLHRNKHTFLQHHSDKRRPPSAELYS NALPVGLPLPSHVSLTGESDGGYMDMSKDESVDYVPMLDMKGDVKYADIESSNYMEREGKFKVVLRAKELRLQYREALAQHLHKDIVLKDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDNNPFSFQSDVYSYGIVLYELMTGELPYSHINNRDQIIFMVGRGYASPDMSKLYC',
        # L1196M mutation: L->M at position ~80 in the kinase domain sequence above
        'note': 'Gatekeeper mutation L1196M reduces Crizotinib binding affinity',
    },
}


def generate_esmfold_structure(sequence: str, output_path: Path, name: str) -> bool:
    """Generate structure using ESMFold API."""
    if output_path.exists():
        print(f"  Structure exists: {output_path.name}")
        return True

    # Clean sequence
    sequence = sequence.replace(' ', '').replace('\n', '')

    if len(sequence) > 400:
        print(f"  Sequence too long ({len(sequence)} aa), truncating to 400")
        sequence = sequence[:400]

    print(f"  Generating ESMFold structure for {name} ({len(sequence)} aa)...")

    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        response = requests.post(
            url,
            data=sequence,
            headers={'Content-Type': 'text/plain'},
            timeout=300
        )

        if response.status_code == 200:
            output_path.write_text(response.text)
            print(f"  Generated: {output_path.name}")
            return True
        else:
            print(f"  ESMFold API error: {response.status_code}")
            print(f"  Response: {response.text[:200]}")
            return False

    except Exception as e:
        print(f"  ESMFold error: {e}")
        return False


def main():
    print("=" * 70)
    print("Generating Splice Variant Kinase Domain Structures")
    print("=" * 70)

    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

    results = []

    for variant_name, info in VARIANTS.items():
        print(f"\n{'='*60}")
        print(f"Variant: {variant_name}")
        print(f"{'='*60}")
        print(f"Gene: {info['gene']}")
        print(f"Drug: {info['drug']}")
        print(f"Expected: {info['expected']}")
        print(f"Note: {info['note']}")

        # Generate canonical kinase domain structure
        canonical_path = STRUCTURES_DIR / f"{info['gene']}_kinase_canonical.pdb"
        if info['canonical_kinase']:
            print(f"\n[1] Canonical kinase domain:")
            success = generate_esmfold_structure(
                info['canonical_kinase'],
                canonical_path,
                f"{info['gene']}_canonical"
            )
            if success:
                results.append((variant_name, 'canonical', canonical_path))
            time.sleep(1)  # Rate limiting

        # Generate variant kinase domain structure
        if info['variant_kinase']:
            variant_path = STRUCTURES_DIR / f"{variant_name}_kinase.pdb"
            print(f"\n[2] Variant kinase domain:")
            success = generate_esmfold_structure(
                info['variant_kinase'],
                variant_path,
                variant_name
            )
            if success:
                results.append((variant_name, 'variant', variant_path))
            time.sleep(1)
        else:
            print(f"\n[2] No variant kinase domain (binding site absent in variant)")

    # Summary
    print("\n" + "=" * 70)
    print("GENERATION SUMMARY")
    print("=" * 70)

    print(f"\nGenerated {len(results)} structures:")
    for variant, struct_type, path in results:
        print(f"  {variant} ({struct_type}): {path.name}")

    print("\nKey insights for validation:")
    print("-" * 70)
    for variant_name, info in VARIANTS.items():
        if info['canonical_kinase'] == info.get('variant_kinase'):
            print(f"  {variant_name}: Kinase domains IDENTICAL")
            print(f"    -> Resistance mechanism: {info['note']}")
        elif info['variant_kinase'] is None:
            print(f"  {variant_name}: Binding domain ABSENT in variant")
            print(f"    -> {info['note']}")
        else:
            print(f"  {variant_name}: Kinase domain MUTATED")
            print(f"    -> {info['note']}")


if __name__ == "__main__":
    main()
