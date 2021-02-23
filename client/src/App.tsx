import * as React from "react";
import { InterpolationPage } from "./InterpolationPage";
import { Navbar, Nav, Form, Button } from "react-bootstrap";
import { DEFAULT_COLLECTION, EActiveTabs, ICollection, IInterpolatedParticle } from "./interfaces";
import { EmbeddingPage } from "./EmbeddingPage";
import { downloadFile } from "./utils";

function App() {
  const [activeTab, setActiveTab] = React.useState<EActiveTabs>(EActiveTabs.EMBEDDING);
  const [collections, setCollections] = React.useState<ICollection[]>([]);
  const [interpolationStructures, setInterpolationStructures] = React.useState<string[]>([
    "NC1CC1c1ccccc1",
    "O=C(CN1C(=O)CSc2ccc(S(=O)(=O)N3CCCCC3)cc21)NCc1cccnc1",
    "COC(=O)[C@H]1C[C@H]2[C@@H]3CCC(=O)[C@@]3(C)CC[C@@H]2[C@@]2(C)CC/C(=NOCCN)C[C@H]12",
  ]);

  React.useEffect(() => {
    collections.forEach((collection) =>
      collection.data.forEach((p, i) => {
        p.index = i;
        p.collection = collection.name;
      })
    );
  }, [collections]);

  const interpolationCollection = React.useMemo(
    () => collections.find((c) => c.name === "Interpolated") as ICollection<IInterpolatedParticle> | undefined,
    [collections]
  );

  const setInterpolationCollection = React.useCallback(
    (collection: ICollection<IInterpolatedParticle>) => {
      setCollections((collections) => [...collections.filter((c) => c.name !== "Interpolated"), collection]);
    },
    [setCollections]
  );

  return (
    <div className="vh-100 d-flex flex-column">
      <Navbar collapseOnSelect expand="lg" bg="dark" variant="dark">
        <Navbar.Brand href="#home">CDDD Explorer</Navbar.Brand>
        <Navbar.Toggle aria-controls="responsive-navbar-nav" />
        <Navbar.Collapse id="responsive-navbar-nav">
          <Nav className="mr-auto">
            <Nav className="mr-auto">
              <Nav.Link
                href="#"
                active={activeTab === EActiveTabs.EMBEDDING}
                onSelect={() => setActiveTab(EActiveTabs.EMBEDDING)}
              >
                Embedding
              </Nav.Link>
              <Nav.Link
                href="#"
                active={activeTab === EActiveTabs.INTERPOLATION}
                onSelect={() => setActiveTab(EActiveTabs.INTERPOLATION)}
              >
                Interpolation
              </Nav.Link>
            </Nav>
          </Nav>
          <Form inline>
            <Button
              variant="outline-primary"
              onClick={() => {
                downloadFile(collections.find((c) => c.name === DEFAULT_COLLECTION)?.data, "export");
              }}
            >
              Export
            </Button>
          </Form>
        </Navbar.Collapse>
      </Navbar>
      <div className="container-fluid mt-2" style={{ flex: 1, overflow: "auto" }}>
        <div className="row" style={{ height: "100%", overflow: "auto", position: "relative" }}>
          {activeTab === EActiveTabs.EMBEDDING ? (
            <EmbeddingPage
              collections={collections}
              setCollections={setCollections}
              setInterpolationStructures={setInterpolationStructures}
              setActiveTab={setActiveTab}
            />
          ) : null}
          {activeTab === EActiveTabs.INTERPOLATION ? (
            <InterpolationPage
              collection={interpolationCollection}
              setCollection={setInterpolationCollection}
              structures={interpolationStructures}
              setStructures={(structures: string[]) => setInterpolationStructures(structures)}
            />
          ) : null}
        </div>
      </div>
    </div>
  );
}

export default App;
