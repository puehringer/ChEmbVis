import {
  ICellRendererFactory,
  ERenderMode,
  ICellRenderer,
  IDataRow,
  IRenderContext,
  IGroupCellRenderer,
  IOrderedGroup,
  renderMissingDOM,
} from "lineupjs";
import { getImageURL, getReducedImages } from "../../utils/api";
import { StructureImageColumn } from "./StructureImageColumn";
import { abortAble } from "lineupengine";

const template =
  '<a target="_blank" rel="noopener" style="background-size: contain; background-position: center; background-repeat: no-repeat;"></a>';

export function svgToImageSrc(svg: string): string {
  return `data:image/svg+xml;base64,${btoa(svg)}`;
}

export function svgToCSSBackground(svg: string): string {
  return `url('${svgToImageSrc(svg)}')`;
}

export class StructureImageRenderer implements ICellRendererFactory {
  readonly title: string = "Chemical Structure";

  canRender(col: StructureImageColumn, mode: ERenderMode): boolean {
    return col instanceof StructureImageColumn && (mode === ERenderMode.CELL || mode === ERenderMode.GROUP);
  }

  create(col: StructureImageColumn): ICellRenderer {
    return {
      template,
      update: (n: HTMLLinkElement, d: IDataRow) => {
        if (!renderMissingDOM(n, col, d)) {
          if(d.v.images?.[0]) {
            n.style.backgroundImage = svgToCSSBackground(d.v.images[0]);
            return;
          }
          const value = col.getValue(d)!;
          // Load aysnc to avoid triggering
          return abortAble(new Promise((resolve) => setTimeout(() => resolve(value), 500))).then((image) => {
            if (typeof image === "symbol") {
              return;
            }
            n.style.backgroundImage = `url('${getImageURL(value, col.getFilter()?.filter, col.getAlign())}')`;
            n.title = value;
            n.href = `https://pubchem.ncbi.nlm.nih.gov/#query=${value}`;
          });
        }
      },
    };
  }

  createGroup(col: StructureImageColumn, context: IRenderContext): IGroupCellRenderer {
    return {
      template,
      update: (n: HTMLImageElement, group: IOrderedGroup) => {
        context.tasks.groupRows(col, group, "StructureImageRendererGroup", (rows) => {
          return abortAble(getReducedImages(Array.from(rows.map((row) => col.getLabel(row))))).then((res: any) => {
            n.style.backgroundImage = res ? svgToCSSBackground(res) : "";
          });
        });
      },
    };
  }
}
